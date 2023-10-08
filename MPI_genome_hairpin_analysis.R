library(Biostrings)
library(AnnotationHub)
library(magrittr)
library(purrr)
library(GenomicRanges)
library(rtracklayer)
library(ensembldb)
library(stringr)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(Rmpi)

# Register the mpi environment for BiocParallel.
# Preserve one core as administrator.
param <- SnowParam(mpi.universe.size() - 1, "MPI")
register(param)

# Get mouse ENSEMBL database version 99.
edb <- query(AnnotationHub(), pattern = c("Mus musculus", "EnsDb", 99))[[1]]
# Get mouse genome sequence.
Mmusculus <- ensembldb:::getGenomeTwoBitFile(edb)

# Get the genes from the ENSEMBL database, then split into GRangesList for later lapply.
mouse_genes_granges <- genes(edb)
mouse_genes_grangeslist <- split(mouse_genes_granges, mouse_genes_granges$gene_id)[names(mouse_genes_granges)]

Mmusculus_seqinfo <- seqinfo(mouse_genes_granges)

# Get the transcript IDs for each mouse gene.
mouse_transcripts <- transcriptsBy(edb, by = "gene")[names(mouse_genes_granges)]
mouse_genes_transcriptIDs <- bplapply(mouse_transcripts, FUN = function(granges) {granges$tx_id})

# Get the coding sequences for each mouse gene.
mouse_genes_cds_grangeslist <- cdsBy(edb, by = "gene")
# Some genes don't have coding sequence.
mouse_genenames_without_cds <- names(mouse_genes_granges)[-which(names(mouse_genes_granges) %in% names(mouse_genes_cds_grangeslist))]
names(mouse_genenames_without_cds) <- mouse_genenames_without_cds
# Assign the missing genes empty CDS GRanges.
mouse_genes_cds_grangeslist <- c(mouse_genes_cds_grangeslist, GRangesList(lapply(mouse_genenames_without_cds, FUN = function(gene_id) {GRanges(seqinfo = seqinfo(mouse_genes_cds_grangeslist))})))[names(mouse_genes_granges)]

# Get the intron, exon, 5'UTR and 3'UTR lists in the same order as the gene granges.
mouse_transcript_exons_grangeslist <- exonsBy(edb, by = "tx")
mouse_genes_exons_grangeslist <- GRangesList(
	bplapply(mouse_genes_transcriptIDs,
		FUN = function(transcriptids) {
			mouse_transcript_exons_grangeslist[transcriptids] %>% unlist()
		}
	)
)[names(mouse_genes_granges)]
rm(mouse_transcript_exons_grangeslist)

mouse_transcript_introns_grangeslist <- intronsByTranscript(edb)
mouse_genes_introns_grangeslist <- GRangesList(
	bplapply(mouse_genes_transcriptIDs,
		FUN = function(transcriptids) {
			mouse_transcript_introns_grangeslist[transcriptids] %>% unlist()
		}
	)
)[names(mouse_genes_granges)]
rm(mouse_transcript_introns_grangeslist)

mouse_transcript_fiveUTR_grangeslist <- fiveUTRsByTranscript(edb)
mouse_genes_fiveUTR_grangeslist <- GRangesList(
	bplapply(mouse_genes_transcriptIDs,
		FUN = function(transcriptids) {
			relevant_ids <- transcriptids[which(transcriptids %in% names(mouse_transcript_fiveUTR_grangeslist))]
			if (length(relevant_ids) == 0) return(GRanges(seqinfo = Mmusculus_seqinfo))
			mouse_transcript_fiveUTR_grangeslist[relevant_ids] %>% unlist()
		}
	)
)[names(mouse_genes_granges)]
rm(mouse_transcript_fiveUTR_grangeslist)

mouse_transcript_threeUTR_grangeslist <- threeUTRsByTranscript(edb)
mouse_genes_threeUTR_grangeslist <- GRangesList(
	bplapply(mouse_genes_transcriptIDs,
		FUN = function(transcriptids) {
			relevant_ids <- transcriptids[which(transcriptids %in% names(mouse_transcript_threeUTR_grangeslist))]
			if (length(relevant_ids) == 0) return(GRanges(seqinfo = Mmusculus_seqinfo))
			mouse_transcript_threeUTR_grangeslist[relevant_ids] %>% unlist()
		}
	)
)[names(mouse_genes_granges)]
rm(mouse_transcript_threeUTR_grangeslist)


#A function for taking a transcript GRanges, finding palindromes, and returning them as genomic GRanges locations
getHairpinsFromMouseGrange <- function(transcript_grange, minimum = 30, mismatch = 3, potential_exons, potential_fiveUTRs, potential_threeUTRs, potential_introns, potential_cds) {

	# Get the palindromes to work with.
	exon_sequences <- getSeq(Mmusculus, transcript_grange)
	transcript_sequence <- unlist(exon_sequences)
	transcript_exon_ends <- cumsum(width(exon_sequences))
	transcript_exon_iranges <- IRanges(end = transcript_exon_ends, width = width(exon_sequences))

	func_palindromes <- Biostrings::findPalindromes(transcript_sequence, min.armlength = minimum, max.looplength = length(transcript_sequence), max.mismatch = mismatch)
	# Return empty list if no palindromes.
	if (length(func_palindromes) == 0) {
		return(GRangesList())
	}

	# Collect all the relevant introns and exons and UTRs for the gene once.
	func_matching_exons <- potential_exons[overlapsAny(potential_exons, transcript_grange, ignore.strand = TRUE)]
	func_matching_fiveUTRs <- potential_fiveUTRs[overlapsAny(potential_fiveUTRs, transcript_grange, ignore.strand = TRUE)]
	func_matching_threeUTRs <- potential_threeUTRs[overlapsAny(potential_threeUTRs, transcript_grange, ignore.strand = TRUE)]
	func_matching_introns <- potential_introns[overlapsAny(potential_introns, transcript_grange, ignore.strand = TRUE)]
	func_matching_cds <- potential_cds[overlapsAny(potential_cds, transcript_grange, ignore.strand = TRUE)]

	# Define the function that gets used repeatedly.
	makeGRangeFromPalindrome <- function(palindrome, original_grange, transcript_exon_iranges, mismatch = 0) {
		stopifnot(length(palindrome) == 1)

		# Get the relationships between the palindrome arms and the transcript.
		left_arm_view <- palindromeLeftArm(palindrome, max.mismatch = mismatch)
		left_arm_exonhits <- subjectHits(findOverlaps(left_arm_view, transcript_exon_iranges))
		left_arm_iranges <- pintersect(findOverlapPairs(left_arm_view, transcript_exon_iranges))
		right_arm_view <- palindromeRightArm(palindrome, max.mismatch = mismatch)
		right_arm_exonhits <- subjectHits(findOverlaps(right_arm_view, transcript_exon_iranges))
		right_arm_iranges <- pintersect(findOverlapPairs(right_arm_view, transcript_exon_iranges))

		if (as.character(strand(original_grange)[1]) == "-") {
			# End(exon_on_genome) + start(exon_on_transcript) - start(palindromearm_on_transcript)
			# 1-based errors cancel each other.
			left_arm_ends <- map_int(left_arm_exonhits, ~ end(original_grange[.x]) + start(transcript_exon_iranges[.x]) - start(intersect(left_arm_view, transcript_exon_iranges[.x])))
			left_arm_widths <- width(left_arm_iranges)
			left_arm_starts <- left_arm_ends - left_arm_widths + 1
			right_arm_ends <- map_int(right_arm_exonhits, ~ end(original_grange[.x]) + start(transcript_exon_iranges[.x]) - start(intersect(right_arm_view, transcript_exon_iranges[.x])))
			right_arm_widths <- width(right_arm_iranges)
			right_arm_starts <- right_arm_ends - right_arm_widths + 1
			# Switch left and right because of strandedness.
			first_arm_starts <- right_arm_starts
			first_arm_ends <- right_arm_ends
			second_arm_starts <- left_arm_starts
			second_arm_ends <- left_arm_ends

		} else {
			# Start(exon_on_genome) - start(exon_on_transcript) + start(palindromearm_on_transcript)
			left_arm_starts <- map_int(left_arm_exonhits, ~ start(original_grange[.x]) - start(transcript_exon_iranges[.x]) + start(intersect(left_arm_view, transcript_exon_iranges[.x])))
			left_arm_widths <- width(left_arm_iranges)
			left_arm_ends <- left_arm_starts + left_arm_widths - 1
			right_arm_starts <- map_int(right_arm_exonhits, ~ start(original_grange[.x]) - start(transcript_exon_iranges[.x]) + start(intersect(right_arm_view, transcript_exon_iranges[.x])))
			right_arm_widths <- width(right_arm_iranges)
			right_arm_ends <- right_arm_starts + right_arm_widths - 1

			first_arm_starts <- left_arm_starts
			first_arm_ends <- left_arm_ends
			second_arm_starts <- right_arm_starts
			second_arm_ends <- right_arm_ends
		}

		# We use "first" and "second" for the arms here so that the left arm of the hairpin is downstream of the right arm on the genome, regardless of the strand of the gene.
		left_arm_granges <- GRanges(
			seqnames = seqnames(original_grange)[1],
			ranges = IRanges(
				start = first_arm_starts,
				end = first_arm_ends
			),
			strand = strand(original_grange)[1],
			seqinfo = seqinfo(original_grange),
			arm = factor("left", levels = c("left", "right"))
		)
		right_arm_granges <- GRanges(
			seqnames = seqnames(original_grange)[1],
			ranges = IRanges(
				start = second_arm_starts,
				end = second_arm_ends
			),
			strand = strand(original_grange)[1],
			seqinfo = seqinfo(original_grange),
			arm = factor("right", levels = c("left", "right"))
		)
		return(sort(c(left_arm_granges, right_arm_granges)))
	}

	checkGRangeGeneRegions <- function(func_grange, potential_exons, potential_fiveUTRs, potential_threeUTRs, potential_introns, potential_cds) {
		func_grange$is_exon <- overlapsAny(func_grange, potential_exons, ignore.strand = TRUE)
		func_grange$is_fiveUTR <- overlapsAny(func_grange, potential_fiveUTRs, ignore.strand = TRUE)
		func_grange$is_threeUTR <- overlapsAny(func_grange, potential_threeUTRs, ignore.strand = TRUE)
		func_grange$is_intron <- overlapsAny(func_grange, potential_introns, ignore.strand = TRUE)
		func_grange$is_cds <- overlapsAny(func_grange, potential_cds, ignore.strand = TRUE)
		return(func_grange)
	}


	return_list <- list()
	#It seems that we can't use lapply() over an XStringViews object, so we must use a for-loop.
	for (i in 1:length(func_palindromes)) {
		return_list[[i]] <- checkGRangeGeneRegions(
			makeGRangeFromPalindrome(func_palindromes[i], transcript_grange, transcript_exon_iranges, mismatch),
			func_matching_exons,
			func_matching_fiveUTRs,
			func_matching_threeUTRs,
			func_matching_introns,
			func_matching_cds
		)
	}
	return(GRangesList(return_list))
}


hairpins_for_genes <- purrr::transpose(
	mcmapply(
		FUN = purrr::safely(getHairpinsFromMouseGrange),
		transcript_grange = mouse_genes_grangeslist,
		potential_exons = mouse_genes_exons_grangeslist,
		potential_fiveUTRs = mouse_genes_fiveUTR_grangeslist,
		potential_threeUTRs = mouse_genes_threeUTR_grangeslist,
		potential_introns = mouse_genes_introns_grangeslist,
		potential_cds = mouse_genes_cds_grangeslist,
		MoreArgs = list(minimum = 30, mismatch = 3),
		SIMPLIFY = FALSE
	)
)
saveRDS(hairpins_for_genes, file = paste0("hairpins_for_mousegenes.safely.", Sys.Date(), ".RDS"))

mpi.quit()
