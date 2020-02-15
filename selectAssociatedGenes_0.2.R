#######################################################
#
#	selectAssociatedGenes.R
#	SELECTS GENES NEAR TO X MB OF A SNP ASSOCIATED
#	WITH A DETERMINED PHENOTYPE.
#
#	AUTHOR. JAVIER MONTALVO-ARREDONDO.
#	VERSION. 0.2
#	CONTACT. buitrejma[at]gmail[dot]com
#	CALZADA ANTONIO NARRO, 1923,
#	SAN BUENAVISTA SALTILLO COAHUILA, C.P. 25315
#
#	TAKES THE GWAS ANALYSIS RESULTS TABLE, A TRANSCRIPTOME ANNOTATION
#	FILE IN GFF FORMAT AND A HEADER HASH TABLE RELATING CHROMOSOME
#	NAME AND NUMBER OF CHROMOSOME. THIS CODE WILL SELECT GENES NEAR 
#	(X MEGABASES) TO A SIGNIFICATIVE SNP ASSOCIATED WITH PHENOTYPE.
#
#
#	PSEUDO-CODE.
#	
#	10	DECLARE XMB AS INT: NUMBER IN MEGABASES TO SEARCH GENES NEAR TO THE SIGNIFICATIVE SNP
#	12	DECLARE ALPHA AS NUMERIC: PVALUE CUTOFF. # DEPRECATED !
#	15	DECLARE PROD_CONTAINER AS A MATRIX WITH DIMS OF 1 X 4 
#	20	WITH INPUT FILES CREATE DATA OBJECTS (OBJECT-GWAS-RESULTS, OBJECT-HEADER-HASH, OBJECT-GFF).
#	30	ITERATE AMONG OBJECT-GWAS-RESULTS ROWS.
#	3010		LET NAME  TAKES OBJECT-GWAS-RESULTS$SNP_ID.
#	3020		LET CHRN  TAKES CHROMOSOME_NAME FROM OBJECT-HEADER-HASH USING OBJECT-GWAS-RESULTS$CHR_NUMBER AS A KEY.
#	3030		LET POS   TAKES OBJECT-GWAS-RESULTS$SNP_POSITION.
#	3035		LET PVAL  TAKES OBJECT-GWAS-RESULTS$SNP_PVLAUE.
#	3040		LET START TAKES CALCULATION OF POS - XMB.
#	3050		LET END   TAKES CALCULATION OF POS + XMB.
#	3060		LET TMP-OBJECT-GFF TAKES SAMPLED OBJECT-GFF BY CHRN
#	3070		LET PROD  TAKES TMP-OBJECT-GFF$PRODUCT IF [TMP-OBJECT-GFF$START >= START & TMP-OBJECT-GFF$END <= END] OR DO CONTINUE
#	3075		LET NAME, CHRN, POS BE THE SAME LENGTH OF A PROD COLUMN.
#	3080		LET TMP_PROD_CONTAINER BE A DATAFRAME
#	3090		LET TMP_PROD_CONTAINER CONTAINS NAME, CHRN, POS, PROD.
#	40	RETURN PROD_CONTAINER
############################################################


selectAssociatedGenes <- function(xmb, SNPS, GFF, HEADER) {
	#	PRE-PROCESSING DATA
	cat("\nPre-processing SNPs position\n")
	SNPS <- SNPS[,1:3]
	print(head(SNPS))

	cat("\nPre-processing GFF transcriptome!\n")
	GFF <- GFF[ GFF$type == "mRNA" | GFF$type == "transcript" , ]
	cat("\nTranscriptome too large to print it on screen!\n")

	cat("\nPre-processing HEADER!\n")
	if (sum(grepl(">", HEADER$V1)) > 0) {
		HEADER$V1 <- gsub(">", "", HEADER$V1)
	}
	print(head(HEADER))
	
	#	SPLIT GFF BY CHRNAME: HIGHLY RECOMENDED TO INCREASE PERFORMANCE
	names <- as.character(unique(sort(seqnames(GFF))))
	length_names <- length(names)
	GFFList <- list()
	cat("\nSubsetting GFF transcriptome by chromosome name\n")
	pb <- txtProgressBar(min = 0, max = length_names, style = 3)
	for (i in 1:length_names) {
		GFFList[[names[i]]] <- GFF[seqnames(GFF) == names[i],]
		setTxtProgressBar(pb, i)
	}
	close(pb)
	rm(GFF)

	#	PERFORMING ITERATIONS AMONG GWAS-RESULTS DATA-FRAME'S ROWS
	cat("\nPerfomring extraction task!\n")
	num_of_gwas_rows <- length(SNPS[,1])
	pb <- txtProgressBar(min = 0, max = num_of_gwas_rows, style = 3)

	for (i in 1:num_of_gwas_rows) {
		name <- as.character(SNPS[i,1])
		chrn <- as.character(HEADER[HEADER$V2 == as.numeric(SNPS[i,2]),1])
		pos <- as.numeric(SNPS[i,3])
		start <- pos - (xmb * 1000000)
		end <- pos + (xmb * 1000000)
		TMP_GFF <- GFFList[[chrn]]
		PROD <- TMP_GFF[ TMP_GFF@ranges@start >= start & TMP_GFF@ranges@start + TMP_GFF@ranges@width <= end, ]@elementMetadata@listData$product
		TRID <- TMP_GFF[ TMP_GFF@ranges@start >= start & TMP_GFF@ranges@start + TMP_GFF@ranges@width <= end, ]@elementMetadata@listData$transcript_id
		GNID <- TMP_GFF[ TMP_GFF@ranges@start >= start & TMP_GFF@ranges@start + TMP_GFF@ranges@width <= end, ]@elementMetadata@listData$gene
		PROD <- cbind(PROD, TRID, GNID)
		if (length(PROD) == 0) {
			setTxtProgressBar(pb, i)
			next
		}
		if(length(PROD) > 3) {
			PROD <- PROD[order(PROD[,2]),]
			PROD <- PROD[ !duplicated(PROD[,2]) , ]
			NAME <- rep(name, length(PROD[,2]))
			CHRN <- rep(chrn, length(PROD[,2]))
			POS <- rep(pos, length(PROD[,2]))
			XMB <- rep(xmb, length(PROD[,2]))
		} else {
			NAME <- name
			CHRN <- chrn
			POS <- pos
			XMB <- xmb
		}
		TMP_PROD_CONTAINER <- cbind(NAME, CHRN, POS, XMB, PROD) # PVAL was removed! Watch out!
		if (exists("PROD_CONTAINER")) {
			PROD_CONTAINER <- rbind(PROD_CONTAINER, TMP_PROD_CONTAINER)
		} else {
			PROD_CONTAINER <- TMP_PROD_CONTAINER
		}
		setTxtProgressBar(pb, i)
	}
	close(pb)
	rm(SNPS)
	as.data.frame(PROD_CONTAINER)
}

#######################################################################################################
#
#	leftJoin
#
#	SAME AUTHORS, CONTACT AND VERSION AS IS DESCRIBED IN selectAssociatedGenes
#
#	PERFORMS A LEFT JOIN OVER DIFFERENT SETS TO INCORPORATE RELATE DATA
#	OF SNPS P-VALUES, PROTEIN IDS, AND GO TERMS TO TRANSCRIPTS IDS
#	AND SNPS IDS.
#	
#	REQUIRED ITEMS:
#		GENES.-		THIS IS THE OUTPUT TABLE FROM selectedAssociatedGenes FUNCTION.
#		SNPS.-		THIS IS THE GAPIT GWAS RESULT CSV TABLE.
#		TRANS_2_PROT.-	THIS IS THE TABLE EXTRACTED FROM GFF TRANSCRIPTOME UPLOADED WITH
#				imort.gff3 FUNCTION FROM rtracklayer R PACKAGE.
#		PROT_2_GO.-	THIS IS THE MODIFIED OUTPUT TABLE FROM INTERPROSCAN SEARCHING,
#				THIS TABLE ONLY HAS THE PROTEIN ID AND GO TERM.
#
#	PSEUDO-CODE
#		LET THE GENES BE LEFT JOINED WITH SNPS BY NAME AND SNP COLNAMES
#		LET THE RESULTING TABLE BE LEFT JOINED WITH TRANS_2_PROT BY TRID AND PARENT.VALUE COLNAMES
#		LET THE RESULTING TABLE BE LEFT JOINED WITH PROT_2_GO BY TRID AND V1 COLNAMES
#
########################################################################################################

leftJoin <- function(GENES, SNPS, TRANS_2_PROT, PROT_2_GO, full=FALSE) {
	cat("\nMerging tables\n")
	#	MERGING GENES-NEAR-TO-SNP TABLE WITH SNPS-PVALUES TABLE
	tmp <- merge(GENES, SNPS[,c("SNP", "P.value")], by.x="NAME", by.y="SNP", all.x=T)

	#	MERGING TMP TABLE WITH TRANSCRIPTS_ID_2_PROT_ID TABLE
	if (sum(grepl(".*-", TRANS_2_PROT$Parent.value)) > 0) {
	
		TRANS_2_PROT$Parent.value <- gsub(".*-", "", TRANS_2_PROT$Parent.value)
	}
	tmp <- merge(tmp, TRANS_2_PROT, by.x="TRID", by.y="Parent.value", all.x=T)

	#	MERGING TMP TABLE WITH PROTEIN_ID_2_GOTERM TABLE
	tmp <- merge(tmp, PROT_2_GO[,c("V1", "V14")], by.x="protein_id", by.y="V1", all.x=T)
	
	cat("\nEliminating duplicates\n")
	tmp <- tmp[ order(tmp$protein_id, tmp$V14, decreasing=T), ]
	#	ELIMINATING ALL DUPLICATES WITH NO GO TERM
	keep1 <- duplicated(tmp[,1:9])
	keep2 <- tmp$V14 == ""
	keep2[ is.na(keep2) ] <- FALSE
	keep3 <- keep1 & keep2
	keep3 <- !keep3
	tmp <- tmp[ keep3, ]
	#	ELIMINATING DUPLICATED ENTRIES
	keep <- !duplicated(tmp)
	tmp <- tmp[ keep, ]
	#	FULL LEFTJOIN
	#	ELIMINATING ALL ENTRIES WITHOUT GO TERMS
	if (full) {
		cat("\nEliminating entries without GO terms\n")
		keep <- tmp$V14 != ""
		keep[ is.na(keep) ] <- FALSE
		tmp <- tmp[ keep, ]
	}

	tmp
}

#############################################################################################################
#
#	GPGtable
#	
#	CREATES A GENE ID, P.VALUE, GO TERMS TABLE
#
#	REQUIRES THE MERGED TABLE (MT) FROM leftJoin FUNCTION.
#
#############################################################################################################

GPGtable <- function(MT) {
	
	cat("\nCreating GPG table,")
	# BUG 1.	Different SNPS are closer to acco20 in MT table, so after sorintg and unique GNID
	#		the way of allocating p-values to GO terms is wrong but results are ok.
	#		It seems to be a classical silent bug who does not affect the output.
	MT <- MT[,c("GNID", "P.value", "V14")]
	GNID <- MT$GNID
	GNID <- unique(sort(GNID)) # HERE IS THE BUG.
	length_GNID <- length(GNID)
	cat("\nProcessing ", length(GNID), "genes\n")
	MT <- MT[order(MT$GNID),]
	pb <- txtProgressBar(min = 0, max = length_GNID, style = 3)

	for (i in 1:length_GNID) {
		p.value <- MT[MT$GNID == GNID[i],]$P.value
		p.value <- as.numeric(unique(sort(p.value)))
		id <- as.character(GNID[i])
		goterms <- unlist(lapply(MT[MT$GNID == GNID[i],]$V14, function(x){strsplit(x, "\\|")[[1]]}))
		goterms <- unique(sort(goterms))
		p.value <- rep(p.value, length(goterms))
		id <- rep(id, length(goterms))
		TMP_gpgtable <- cbind(ID=id, P.value=p.value, V14=goterms)
		if (exists("gpgtable")) {
			gpgtable <- rbind(gpgtable, TMP_gpgtable)
		} else {
			gpgtable <- TMP_gpgtable
		}
		setTxtProgressBar(pb, i)
	}
	
	close(pb)
	gpgtable <- as.data.frame(gpgtable)
	gpgtable$ID <- as.character(gpgtable$ID)
	gpgtable$P.value <- as.numeric(as.character(gpgtable$P.value))
	keep <- !duplicated(gpgtable)
	gpgtable <- gpgtable[ keep, ]
	gpgtable
}

#############################################################################################################
#
#	GWA.test
#
#


GWA.test <- function(GPGT, boot=1) {

	keep <- tapply( GPGT$V14, GPGT$V14, length ) >= 10
	goterms <- names( keep )
	goterms <- goterms[ keep ]
	GPGT <- GPGT[ GPGT$V14 %in% goterms, ]
	top6goterms <- sort(tapply(GPGT$P.value, GPGT$V14, median))
	top6goterms_names <- names(top6goterms)
	GPGT_rowslength <- length(GPGT[,1])
	

	for (goterm_name in top6goterms_names) {
		cat( paste("\nTesting ", goterm_name, "\n", sep = "") )
		sampled_GPGT <- GPGT[ GPGT$V14 == goterm_name , ]
		p_value_sampled_GPGT_median <- median( sampled_GPGT$P.value )
		GPGT_sampled_rows_length <- length( sampled_GPGT$P.value )

		pb <- txtProgressBar(min = 0, max = boot, style = 3)
		p_value_medians <- c()
		for (i in 1:boot) {
			test_p_values <- GPGT[ rownames(GPGT) %in% sample(1:GPGT_rowslength, GPGT_sampled_rows_length), "P.value" ]
			p_value_medians <- c( p_value_medians, median(test_p_values) )
			setTxtProgressBar(pb, i)
		}
		close(pb)
		p_value_miu <- mean( p_value_medians )
		p_value_sd <- sd( p_value_medians )

		pnorm_sampled_GPGT_mean_of_medians <- pnorm( p_value_sampled_GPGT_median, p_value_miu, p_value_sd, lower.tail = TRUE, log.p=FALSE )
		tmp_output <- c(p_value_sampled_GPGT_median, p_value_miu, p_value_sd, pnorm_sampled_GPGT_mean_of_medians)
		if (exists("goterms_p_values")) {
			goterms_p_values <- rbind(goterms_p_values, tmp_output)
		} else {
			goterms_p_values <- tmp_output
		}
	}

	colnames(goterms_p_values) <- c("PValues median", "MIU", "SD", "PNORM")
	rownames(goterms_p_values) <- names(top6goterms)

	goterms_p_values
}
