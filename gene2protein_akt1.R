##############################################################################
# Gene2Protein Assignment Solution
# Target Gene: AKT1 (ENSG00000142208)
# Environment: VSCode + Conda R Environment
##############################################################################

############################
# Load Required Packages
############################

suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(IRanges))

cat("\nInitializing AnnotationHub...\n")

############################
# Connect to AnnotationHub
############################

# Increase timeout for slow internet
options(timeout = 300)
options(download.file.method = "libcurl")

cat("Connecting to AnnotationHub...\n")

ah <- tryCatch({

    # Try online connection first
    AnnotationHub(localHub = FALSE)

}, error = function(e) {

    cat("\n[WARNING] Network connection failed.\n")
    cat("Details:\n")
    cat(e$message, "\n")

    cat("Attempting local cache access...\n")

    # Fallback to local cache
    AnnotationHub(localHub = TRUE)
})

cat("Connected successfully.\n")

############################
# Load Ensembl Database
############################

cat("\nFetching EnsDb database...\n")

edb <- tryCatch({

    # Professor specified database
    ah[["AH119325"]]

}, error = function(e) {

    cat("\n[WARNING] AH119325 unavailable.\n")
    cat("Searching for latest Homo sapiens EnsDb...\n")

    human_ensdbs <- query(ah, c("EnsDb", "Homo sapiens"))

    if (length(names(human_ensdbs)) == 0) {

        stop(
            "FATAL ERROR: No EnsDb database found.\n",
            "Run script once with active internet."
        )
    }

    # Use latest database
    fallback_id <- tail(names(human_ensdbs), 1)

    cat(sprintf("Using fallback database: %s\n", fallback_id))

    ah[[fallback_id]]
})

############################
# Define Target Gene
############################

target_gene_id   <- "ENSG00000142208"
target_gene_name <- "AKT1"

cat(sprintf(
    "\nStarting analysis for %s (%s)\n",
    target_gene_name,
    target_gene_id
))

##############################################################################
# 1. FETCH ALL TRANSCRIPTS
##############################################################################

cat("\nFetching transcripts...\n")

txs <- transcripts(
    edb,
    filter = ~ gene_id == target_gene_id
)

# Create dataframe
tx_df <- data.frame(
    Transcript_ID = txs$tx_id
)

# Save transcripts
write.table(
    tx_df,
    file = "AKT1_transcripts.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat(sprintf(
    "Saved %d transcripts to AKT1_transcripts.tsv\n",
    nrow(tx_df)
))

##############################################################################
# 2. FETCH ALL PROTEINS
##############################################################################

cat("\nFetching proteins...\n")

prts <- proteins(
    edb,
    filter = ~ gene_id == target_gene_id
)

# Remove NA or empty IDs
protein_ids <- prts$protein_id[
    !is.na(prts$protein_id) &
    prts$protein_id != ""
]

# Keep unique proteins
protein_ids <- unique(protein_ids)

# Create dataframe
prt_df <- data.frame(
    Protein_ID = protein_ids
)

# Save proteins
write.table(
    prt_df,
    file = "AKT1_proteins.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat(sprintf(
    "Saved %d proteins to AKT1_proteins.tsv\n",
    nrow(prt_df)
))

##############################################################################
# 3. PROTEIN TO GENOME MAPPING
##############################################################################

cat("\nMapping proteins to genome coordinates...\n")

mapping_results <- data.frame()

# Process each protein
for (pid in protein_ids) {

    cat(sprintf("Processing protein: %s\n", pid))

    # Fetch protein data
    protein_data <- proteins(
        edb,
        filter = ~ protein_id == pid
    )

    protein_df <- as.data.frame(protein_data)

    # Skip empty result
    if (nrow(protein_df) == 0) {

        cat(" -> No protein data found\n")
        next
    }

    # Check for protein_length column
    if (!"protein_length" %in% colnames(protein_df)) {

        cat(" -> protein_length column missing\n")
        next
    }

    protein_length <- protein_df$protein_length[1]

    # Validate protein length
    if (
        length(protein_length) == 0 ||
        is.na(protein_length) ||
        protein_length <= 0
    ) {

        cat(" -> Invalid protein length\n")
        next
    }

    # Create IRanges object
    ir <- IRanges(
        start = 1,
        end = protein_length
    )

    # Run protein to genome mapping
    mapped <- tryCatch({

        proteinToGenome(
            ir,
            edb,
            protein_id = pid
        )

    }, error = function(e) {

        cat(sprintf(
            " -> Mapping failed: %s\n",
            e$message
        ))

        return(NULL)
    })

    # Skip failed mappings
    if (is.null(mapped) || length(mapped) == 0) {
        next
    }

    # Convert to dataframe
    tmp_df <- as.data.frame(mapped)

    # Skip empty dataframe
    if (nrow(tmp_df) == 0) {
        next
    }

    # Add protein ID
    tmp_df$Protein_ID <- pid

    # Append results
    mapping_results <- rbind(
        mapping_results,
        tmp_df
    )
}

############################
# Save Mapping Results
############################

if (nrow(mapping_results) > 0) {

    clean_mapping <- data.frame(
        Protein_ID    = mapping_results$Protein_ID,
        Chromosome    = mapping_results$seqnames,
        Genomic_Start = mapping_results$start,
        Genomic_End   = mapping_results$end,
        Strand        = mapping_results$strand
    )

    write.table(
        clean_mapping,
        file = "AKT1_protein_to_genome_map.tsv",
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("\nSaved mapping file: AKT1_protein_to_genome_map.tsv\n")

} else {

    cat("\nWARNING: No mappings generated.\n")
}

##############################################################################
# FINAL MESSAGE
##############################################################################

cat("\n============================================================\n")
cat("Analysis Complete!\n")
cat("Generated Files:\n")
cat("1. AKT1_transcripts.tsv\n")
cat("2. AKT1_proteins.tsv\n")
cat("3. AKT1_protein_to_genome_map.tsv\n")
cat("============================================================\n")

cat("\nNext Steps:\n")
cat("1. Open IGV locally\n")
cat("2. Set genome reference to GRCh38 / hg38\n")
cat("3. Search for AKT1\n")
cat("4. Compare coordinates with generated mapping file\n\n")