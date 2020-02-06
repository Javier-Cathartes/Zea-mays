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
#	TO A SIGNIFICATIVE SNP ASSOCIATED WITH PHENOTYPE.
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
