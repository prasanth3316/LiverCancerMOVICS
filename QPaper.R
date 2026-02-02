if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

BiocManager::install("DESeq2",force=TRUE)
BiocManager::install("ConsensusClusterPlus",force=TRUE)
Sys.setenv('R_MAX_VSIZE' = '32G')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# --- 2. Install TCGAbiolinks, if you don't have it ---
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

# --- 1. Install dplyr (if you don't have it) ---
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

# --- 2. Load the library ---
# This is the line that fixes your error
library(dplyr)

# --- 3. Load the library ---
# This is the line that fixes your error
library(TCGAbiolinks)
library(DESeq2)
library(pheatmap)
library(ConsensusClusterPlus)
library(RColorBrewer) # For heatmap colors

# 1.1. Define the query for *only* primary tumor samples
query_exp <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - FPKM",
  sample.type = "Primary Tumor"  # This filters *before* you download
)

# 1.2. Download the files
# This will be a smaller download because it's only tumor samples
try(GDCdownload(query_exp, method = "api", files.per.chunk = 50))

# 1.3. Prepare the data
# This creates a SummarizedExperiment object
data_exp_se <- GDCprepare(query_exp)


# 2.1. Extract the raw count matrix
# We use "unstranded" for raw counts from the STAR-Counts workflow
# Rows = all gene types (Ensembl IDs), Columns = full sample barcodes
all_gene_counts_matrix <- assay(data_exp_se, "unstranded")

# 2.2. Extract the gene annotation data
# This tells us which gene is which (e.g., protein_coding, lncRNA)
gene_annotation_df <- as.data.frame(rowData(data_exp_se))

# 2.3. Get the list of *only* protein-coding (mRNA) genes
mrna_gene_ids <- gene_annotation_df %>%
  filter(gene_type == "protein_coding") %>%
  pull(gene_id) # 'gene_id' is the Ensembl ID

# 2.4. Filter your matrix for *only* those genes
# This is your "rows as genes" step
mrna_expression_matrix <- all_gene_counts_matrix[rownames(all_gene_counts_matrix) %in% mrna_gene_ids, ]

# 2.5. Rename columns to 12-character patient ID
# This is your "columns as patient id" step
# Because we queried for 'Primary Tumor', there are no duplicates
colnames(mrna_expression_matrix) <- substr(colnames(mrna_expression_matrix), 1, 12)

# 2.6. Done! View your final matrix
print("--- Gene Expression Matrix (mRNA, Tumor Only) ---")
print(paste("Dimensions:", dim(mrna_expression_matrix)[1], "genes,", dim(mrna_expression_matrix)[2], "patients"))
print(mrna_expression_matrix[1:5, 1:5])


# --- 1. Start with your 'mrna_expression_matrix' ---
# (This object must be in your R environment)
print(paste("Original data dimensions:", dim(mrna_expression_matrix)[1], "genes,", dim(mrna_expression_matrix)[2], "patients"))


# --- 2. Normalize and Filter (as requested) ---

# 2.1. Create a DESeq2 object
# We need a simple 'colData' file
col_data <- data.frame(row.names = colnames(mrna_expression_matrix), 
                       patient_id = colnames(mrna_expression_matrix))

# Create the object (using a 'dummy' design: ~ 1)
# We must round the counts to integers for DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(mrna_expression_matrix), 
                              colData = col_data,
                              design = ~ 1)

# 2.2. Apply Variance Stabilizing Transformation (vst)
# This normalizes the data, making it suitable for clustering
print("Running VST normalization... (This may take a minute)")
vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

# 2.3. Select Top Variable Genes (your "select top genes variables" step)
# We will select the top 1500 most variable genes
print("Selecting top 1500 variable genes...")
rv <- rowVars(vst_matrix)
select <- order(rv, decreasing = TRUE)[1:1500]
matrix_for_clustering <- vst_matrix[select, ]

print(paste("Matrix for clustering dimensions:", dim(matrix_for_clustering)[1], "genes,", dim(matrix_for_clustering)[2], "patients"))

results_dir <- "consensus_clustering_mrna"
if (!dir.exists(results_dir)) dir.create(results_dir)

# 2.2. Run the Consensus Clustering
print("Running Consensus Clustering for k=2 to k=6...")
cc_results <- ConsensusClusterPlus(
  matrix_for_clustering,  # Use the normalized, filtered data
  maxK = 6,               # Test k=2 through k=6
  reps = 100,             # Number of bootstraps (use 100 for a test, 1000 for a paper)
  pItem = 0.8,            # Use 80% of samples in each bootstrap
  pFeature = 1,           # Use 100% of the (already filtered) genes
  clusterAlg = "pam",     # "pam" (k-medoids) is robust
  distance = "euclidean",
  title = results_dir,    # Output directory name
  plot = "pdf",           # Save plots as PDF
  seed = 12345            # For reproducible results
)

print(paste("Consensus clustering complete. Plots saved to:", results_dir))

# --- Let's say you decided k=3 was the best ---
optimal_k <- 3

# Get cluster assignments for k=3
patient_clusters <- cc_results[[optimal_k]]$consensusClass

# This gives you a named vector:
# TCGA-BC-A10Q TCGA-BC-A10R TCGA-BC-A10S ...
#            1            3            1 ...

print("Patient cluster assignments (for k=3):")
print(head(patient_clusters))

# You can now save this as a dataframe
cluster_df <- data.frame(
  patient_id = names(patient_clusters),
  cluster = patient_clusters
)

print(head(cluster_df))

# --- 1.1. Start with your 'mrna_expression_matrix' ---
print(paste("Original data dimensions:", dim(mrna_expression_matrix)[1], "genes,", dim(mrna_expression_matrix)[2], "patients"))

# --- 1.2. Normalize the data (using VST) ---
# Create a 'colData' file for DESeq2
col_data <- data.frame(row.names = colnames(mrna_expression_matrix), 
                       patient_id = colnames(mrna_expression_matrix))

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(mrna_expression_matrix),
                              colData = col_data,
                              design = ~ 1)

# Run VST
print("Running VST normalization...")
vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

# --- 1.3. Filter for Top 500 Variable Genes ---
print("Selecting top 500 variable genes...")
rv <- rowVars(vst_matrix)
select <- order(rv, decreasing = TRUE)[1:500] # Get top 500
matrix_top500 <- vst_matrix[select, ]

print(paste("Matrix for heatmap/clustering dimensions:", dim(matrix_top500)[1], "genes,", dim(matrix_top500)[2], "patients"))

# pheatmap will scale the rows (genes) so you can see the relative
# high/low expression, and it will cluster both rows and columns.
print("Generating heatmap...")

pheatmap(
  matrix_top500,
  show_colnames = FALSE, # Too many patients to show
  show_rownames = FALSE, # Too many genes to show
  main = "Heatmap of Top 500 Variable Genes (mRNA)",
  scale = "row",         # This z-scores the genes (rows)
  cluster_rows = TRUE,   # Cluster the genes
  cluster_cols = TRUE    # Cluster the patients
)


# 3.1. Transpose the matrix
# Rows = patients, Columns = genes
data_for_scoring <- t(matrix_top500)

# 3.2. Scale the data (features/genes)
# This is recommended for distance-based clustering
data_for_scoring_scaled <- scale(data_for_scoring)

print("Calculating internal clustering scores (k=2 to 6)...")

# 3.3. Create Silhouette Score Plot
# This plot shows the "Average Silhouette Width"
# A higher score is better. Look for the peak.
# 1. Install the package (if you haven't already)
install.packages("factoextra")

# 2. Load the package into your R session (this is the crucial step)
library(factoextra)

# 3. Now, re-run your plotting code. It should work.
plot_sil <- fviz_nbclust(
  data_for_scoring_scaled, 
  FUNcluster = pam,      
  method = "silhouette", 
  k.max = 6              
) + 
  labs(title = "Optimal Clusters: Average Silhouette Score") +
  theme_minimal()

print(plot_sil)

# --- 1.1. Start with your 'mrna_expression_matrix' ---
print(paste("Original data dimensions:", dim(mrna_expression_matrix)[1], "genes,", dim(mrna_expression_matrix)[2], "patients"))

# --- 1.2. Normalize the data (using VST) ---
# Create a 'colData' file for DESeq2
col_data <- data.frame(row.names = colnames(mrna_expression_matrix), 
                       patient_id = colnames(mrna_expression_matrix))

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(mrna_expression_matrix),
                              colData = col_data,
                              design = ~ 1)

# Run VST
print("Running VST normalization...")
vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)
find_gaps <- function(annotation_df, variable) {
  gaps <- which(annotation_df[[variable]][-1] != annotation_df[[variable]][-nrow(annotation_df)])
  return(gaps)
}
# --- 1.3. Filter for Top 500 Variable Genes ---
print("Selecting top 500 variable genes...")
rv <- rowVars(vst_matrix)
select <- order(rv, decreasing = TRUE)[1:500] # Get top 500
matrix_top500 <- vst_matrix[select, ]

print(paste("Matrix for heatmap/clustering dimensions:", dim(matrix_top500)[1], "genes,", dim(matrix_top500)[2], "patients"))

# --- 1.4. Prepare data FOR CLUSTERING ---
# We must transpose (patients as rows) and scale (genes as columns)
# This is what the pam() algorithm expects
data_for_clustering <- t(matrix_top500)
data_for_clustering_scaled <- scale(data_for_clustering)
# Define your "red and blue" color scheme
my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Loop from k=2 to k=6
for (k_value in 2:6) {
  
  print(paste("--- Generating heatmap for k =", k_value, "---"))
  
  # 1. Cluster the patients into k groups
  pam_fit <- pam(data_for_clustering_scaled, k = k_value)
  
  # 2. Create the column annotation dataframe
  # This tells pheatmap which cluster each patient belongs to
  annotation_col <- data.frame(
    Cluster = as.factor(pam_fit$clustering)
  )
  rownames(annotation_col) <- colnames(matrix_top500)
  
  # 3. Sort the annotation by cluster number
  # This is the key to creating the "block" heatmap
  annotation_col_sorted <- annotation_col %>%
    arrange(Cluster)
  
  # 4. Get the patient IDs in the new sorted order
  patient_order <- rownames(annotation_col_sorted)
  
  # 5. Generate the heatmap
  pheatmap(
    matrix_top500[, patient_order],   # Sort the matrix columns by cluster
    color = my_colors,                # Use red/blue colors
    scale = "row",                    # Normalize genes (rows), as requested
    cluster_rows = TRUE,              # Cluster the genes (rows)
    cluster_cols = FALSE,             # DO NOT cluster columns (we already sorted them)
    show_colnames = FALSE,            # Hide patient IDs
    show_rownames = FALSE,            # Hide gene IDs
    main = paste("Heatmap of Top 500 Genes (k =", k_value, "Clusters)"),
    annotation_col = annotation_col_sorted,  # Add the cluster annotation bar
    gaps_col = find_gaps(annotation_col_sorted, "Cluster") # Add vertical lines between clusters
  )
}

# Helper function to find where to draw lines between clusters
find_gaps <- function(annotation_df, variable) {
  gaps <- which(annotation_df[[variable]][-1] != annotation_df[[variable]][-nrow(annotation_df)])
  return(gaps)
}
  # --- Start with your scaled, transposed data ---
  # (data_for_scoring_scaled)
  
  print("Calculating Gap Statistic... (This may take a few minutes)")
  
  # 1.1. Create the Gap Statistic plot
  # We use the same PAM clustering algorithm for consistency
  plot_gap <- fviz_nbclust(
    data_for_scoring_scaled,
    FUNcluster = pam,      # Use PAM (K-medoids) algorithm
    method = "gap_stat",   # The score to use
    k.max = 6,             # Test k=1 to k=6
    nboot = 500,            # Number of bootstraps. Increase to 500 for better accuracy.
    verbose = FALSE        # Hide the progress bar
  ) +
    labs(title = "Optimal Clusters: Gap Statistic") +
    theme_minimal()
  
  # 1.2. Show the plot
  print(plot_gap)

  cluster_1_ids <- cluster_df$patient_id[cluster_df$cluster == 1]
  
  # --- Get IDs for Cluster 2 ---
  cluster_2_ids <- cluster_df$patient_id[cluster_df$cluster == 2]
  
  # --- Get IDs for Cluster 3 ---
  cluster_3_ids <- cluster_df$patient_id[cluster_df$cluster == 3]
  
  
  # --- Print the list of patient IDs for Cluster 1 ---
  print("--- Cluster 1 Patient IDs ---")
  cat(cluster_3_ids, sep = "\n")
  
  # 1. Get lncRNA gene IDs
  lncrna_gene_ids <- gene_annotation_df %>%
    filter(gene_type == "lncRNA") %>%
    pull(gene_id)
  
  # 2. Filter your matrix for *only* those lncRNA genes
  lncrna_expression_matrix <- all_gene_counts_matrix[rownames(all_gene_counts_matrix) %in% lncrna_gene_ids, ]
  
  # 3. Rename columns to 12-character patient ID
  colnames(lncrna_expression_matrix) <- substr(colnames(lncrna_expression_matrix), 1, 12)
  
  print("--- lncRNA Expression Matrix (Tumor Only) ---")
  print(lncrna_expression_matrix[1:5, 1:5])
  
  
  print("Querying for Somatic Mutation data...")
  query_mut <- GDCquery(
    project = "TCGA-LIHC",
    data.category = "Simple Nucleotide Variation", 
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking" 
  )
  
  try(GDCdownload(query_mut))
  data_maf_path <- GDCprepare(query_mut)
  
  # Load into maftools
  lihc_maf_object_all <- read.maf(maf = data_maf_path)
  
  # Filter for *only* Primary Tumor barcodes ("TP")
  all_maf_barcodes <- lihc_maf_object_all@clinical.data$Tumor_Sample_Barcode
  tumor_maf_barcodes <- TCGAquery_SampleTypes(barcode = all_maf_barcodes, typesample = "TP")
  maf_tumor_only <- subsetMaf(maf = lihc_maf_object_all, tsb = tumor_maf_barcodes)
  
  # This is the 'maf' object we will use for MOVICS
  print("--- Mutation MAF Object (Tumor Only) ---")
  print(maf_tumor_only)
  
  # 1. Extract the main data (the "long" list) from your MAF object
  maf_data_long <- maf_tumor_only@data
  
  # 2. Create the binary matrix using dplyr/tidyr
  print("Converting MAF object to binary mutation matrix...")
  
  binary_mutation_matrix <- maf_data_long %>%
    
    # Select only the gene name and patient barcode
    select(gene = Hugo_Symbol, patient_barcode = Tumor_Sample_Barcode) %>%
    
    # Shorten the barcode to the 12-character patient ID
    mutate(patient_id = substr(patient_barcode, 1, 12)) %>%
    
    # Select only the columns we need
    select(gene, patient_id) %>%
    
    # Get only unique gene-patient pairs
    # This makes it binary (1 event is the same as 5 events)
    distinct() %>%
    
    # Add a "value" column with 1 for all
    mutate(value = 1) %>%
    
    # Pivot the data "wide"
    # Genes will be rows, patient_id will be columns
    pivot_wider(
      names_from = patient_id,
      values_from = value,
      values_fill = 0  # Patients with no mutation in a gene get 0
    ) %>%
    
    # Convert to a standard data.frame
    as.data.frame() %>%
    
    # Set the gene names as the row names
    column_to_rownames(var = "gene") %>%
    
    # Convert to a matrix (as MOVICS prefers)
    as.matrix()
  
  
  # 3. View your new matrix
  print("--- Binary Mutation Matrix (Genes x Patients) ---")
  print(paste("Dimensions:", dim(binary_mutation_matrix)[1], "genes,", dim(binary_mutation_matrix)[2], "patients"))
  print(binary_mutation_matrix[1:5, 1:5])
  
  # 4. IMPORTANT: Keep your 'maf_tumor_only' object!
  # You will still need the original 'maf_tumor_only' object for maftools functions.
  # The 'binary_mutation_matrix' is what you will feed into MOVICS.
  # --- Load Required Libraries ---
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  
  # --- 1. Query for DNA Methylation (450k) ---
  print("Querying for DNA Methylation data...")
  
  query_meth <- GDCquery(
    project = "TCGA-LIHC",
    data.category = "DNA Methylation",
    platform = "Illumina Human Methylation 450", # Standard TCGA platform
    data.type = "Methylation Beta Value",     # This is the processed data (0 to 1)
    sample.type = "Primary Tumor"             # Filter BEFORE downloading
  )
  
  # --- 2. Download the Data ---
  # This will download only the tumor files (approx. 377 files)
  try(GDCdownload(query_meth, method = "api", files.per.chunk = 20))
  
  # --- 3. Prepare the Data ---
  # This combines all files into one SummarizedExperiment object
  print("Preparing methylation data...")
  data_meth_se <- GDCprepare(query_meth)
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("sesame")
  # --- 4. Generate the Matrix ---
  # Extract the Beta Value matrix from the object
  # Rows = CpG Probes, Columns = Full Sample Barcodes
  methylation_matrix <- assay(data_meth_se)
  
  # Rename columns to 12-character patient ID
  colnames(methylation_matrix) <- substr(colnames(methylation_matrix), 1, 12)
  
  # --- 5. View Your Final Matrix ---
  print("--- Methylation Matrix (Tumor Only) ---")
  print(paste("Dimensions:", dim(methylation_matrix)[1], "probes,", dim(methylation_matrix)[2], "patients"))
  print(methylation_matrix[1:5, 1:5])
  
  # --- 1. Calculate variance for every probe (row) ---
  # We must use na.rm = TRUE to ignore the NAs
  print("Calculating variance for all probes...")
  probe_variances <- apply(methylation_matrix, 1, var, na.rm = TRUE)
  
  # --- 2. Sort the variances (highest first) ---
  sorted_variances <- sort(probe_variances, decreasing = TRUE)
  
  # --- 3. Get the names of the top 1000 probes ---
  top_1000_probes <- names(sorted_variances[1:1000])
  
  # --- 4. Filter your matrix for *only* these 1000 probes ---
  methylation_top1000_matrix <- methylation_matrix[top_1000_probes, ]
  
  print(paste("New matrix dimensions:", dim(methylation_top1000_matrix)[1], "probes,", dim(methylation_top1000_matrix)[2], "patients"))
  print(methylation_top1000_matrix[1:5, 1:5])
 
  # --- 1. Install/Load the imputation package ---
  if (!requireNamespace("impute", quietly = TRUE)) {
    BiocManager::install("impute")
  }
  library(impute)
  
  # --- 2. Run k-NN Imputation ---
  # This fills any remaining NAs in your top 1000 matrix
  print("Running k-NN imputation on top 1000 probes...")
  imputed_data <- impute.knn(methylation_top1000_matrix, k = 10)
  
  # The final, clean matrix is in the 'data' slot
  methylation_matrix_clean <- imputed_data$data
  
  # --- 3. Check your final matrix ---
  print("--- Clean & Imputed Top 1000 Methylation Matrix ---")
  print(methylation_matrix_clean[1:5, 1:5])
  print(paste("Any NAs left?", any(is.na(methylation_matrix_clean))))
  # --- Load Required Libraries ---
  
  
  
  ### =========================================================================
  ### Step 1: Install/Load Packages
  ### =========================================================================
  # We'll use 'data.table' for fast file reading and 'dplyr' for easy data handling.
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  
  library(data.table)
  library(dplyr)
  
  ### =========================================================================
  ### Step 2: Define File URLs (from your prompt)
  ### =========================================================================
  
  # 1. The main data file (genes x samples)
  data_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.gene-level_ascat2.tsv.gz"
  data_file <- "TCGA-LIHC.gene-level_ascat2.tsv.gz"
  
  # 2. The gene mapping file (Ensembl ID -> Gene Symbol)
  # Note: This uses gencode.v36, as it's the one Xena provides for this dataset.
  map_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap"
  map_file <- "gencode.v36.probemap.tsv"
  
  ### =========================================================================
  ### Step 3: Download and Load Files
  ### =========================================================================
  
  # Download the data and map files
  print("Downloading data file...")
  download.file(data_url, data_file)
  
  print("Downloading gene map file...")
  download.file(map_url, map_file)
  
  # Load the files using data.table::fread (it's fast and handles .gz)
  print("Loading data...")
  cna_data_raw <- fread(data_file)
  gene_map <- fread(map_file)
  
  print("Files loaded successfully.")
  
  ### =========================================================================
  ### Step 4: Map Ensembl IDs to Gene Symbols
  ### =========================================================================
  
  # 1. Prepare the gene map (we only need the Ensembl 'id' and 'gene' symbol)
  gene_map_simple <- gene_map %>% 
    select(id, gene) %>% 
    distinct() # Keep only unique mappings
  
  # 2. Rename the first column of our data (which is Ensembl IDs) to 'id' to match the map
  setnames(cna_data_raw, old = 1, new = "id")
  
  # 3. Merge the gene map with the data
  # This matches rows in 'cna_data_raw' to their 'gene' symbol from 'gene_map_simple'
  cna_merged <- merge(gene_map_simple, cna_data_raw, by = "id")
  
  # 4. Handle duplicate Gene Symbols
  # Sometimes multiple Ensembl IDs map to the same gene symbol.
  # We will take the *mean* copy number for these cases, a standard practice.
  cna_aggregated <- cna_merged %>%
    select(-id) %>%                           # Remove the Ensembl ID column
    group_by(gene) %>%                        # Group by Gene Symbol
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% # Calculate mean for all sample columns
    filter(!is.na(gene) & gene != "")         # Remove any blank gene names
  
  print("Gene IDs mapped and aggregated.")
  
  ### =========================================================================
  ### Step 5: Create the Final Matrix for MOVICS
  ### =========================================================================
  
  # Convert the aggregated data.frame into the (genes x samples) matrix MOVICS needs
  
  # 1. Convert from tibble to a standard data.frame
  cna_df <- as.data.frame(cna_aggregated)
  
  # 2. Set the 'gene' column as the row names
  rownames(cna_df) <- cna_df$gene
  
  # 3. Remove the 'gene' column, leaving a pure numeric matrix
  cna_matrix <- as.matrix(cna_df[, -1])
  
  # --- Done! ---
  colnames(cna_matrix) <- substr(colnames(cna_matrix), 1, 12)
  
  # --- 5. Check the final matrix again ---
  print("CNA Matrix is ready for MOVICS.")
  print("Dimensions:")
  print(dim(cna_matrix))
  print("Example data (with 12-character IDs):")
  print(cna_matrix[1:5, 1:5])
  
  ### =========================================================================
  ### Next Step
  ### =========================================================================
  #
  # You can now add this 'cna_matrix' to your 'mo.data' list:
  #
  # mo.data <- list(
  #   mRNA        = mrna_matrix,
  #   lncRNA      = lncrna_matrix,
  #   methylation = methylation_matrix,
  #   cna         = cna_matrix,  # <-- Your new Xena matrix
  #   mutation    = mutation_matrix
  # )
  
  
  
  
  ids_mrna <- colnames(mrna_expression_matrix)
  ids_lncrna <- colnames(lncrna_expression_matrix)
  ids_mut <- colnames(binary_mutation_matrix)
  ids_meth <- colnames(methylation_matrix_clean)
  ids_cna <- colnames(cna_matrix)
  
  # --- 2. Check the number of patients in each dataset ---
  # This is a crucial debugging step!
  print("--- Patients per Omics Dataset (Before Intersection) ---")
  print(paste("mRNA patients:", length(ids_mrna)))
  print(paste("lncRNA patients:", length(ids_lncrna)))
  print(paste("Mutation patients:", length(ids_mut)))
  print(paste("Methylation patients:", length(ids_meth)))
  print(paste("CNA patients:", length(ids_cna)))
  
  # --- 3. Create a list of all your ID vectors ---
  list_of_ids <- list(
    mRNA = ids_mrna,
    lncRNA = ids_lncrna,
    Mutation = ids_mut,
    Methylation = ids_meth,
    cna=ids_cna
  )
  
  # --- 4. Find the intersection (the common patients) ---
  common_patients <- Reduce(intersect, list_of_ids)
  
  # --- 5. See your final patient count ---
  print("-----------------------------------------------------")
  print(paste("Total common patients with all 5 data types:", length(common_patients)))
  print("First 10 common patient IDs:")
  print(head(common_patients, 10))
  #
  
  
  
  
  # --- 1. Define a VST helper function ---
  # This function takes a raw count matrix and returns a normalized matrix
  run_vst <- function(count_matrix) {
    print(paste("Running VST on matrix with", ncol(count_matrix), "samples..."))
    
    # Create the simple 'colData' DESeq2 needs
    col_data <- data.frame(row.names = colnames(count_matrix),
                           patient_id = colnames(count_matrix))
    
    # Create the DESeqDataSet (rounding counts to integers)
    # You must have DESeq2 loaded: library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                  colData = col_data,
                                  design = ~ 1)
    
    # Run VST
    vst_data <- vst(dds, blind = TRUE)
    
    # Return the normalized matrix
    return(assay(vst_data))
  }
  
  # --- 2. Filter your count matrices for ONLY the 354 common patients ---
  mrna_counts_filtered <- mrna_expression_matrix[, common_patients]
  lncrna_counts_filtered <- lncrna_expression_matrix[, common_patients]
  
  # --- 3. Run VST on your filtered matrices ---
  print("Normalizing mRNA data...")
  mrna_norm_matrix <- run_vst(mrna_counts_filtered)
  
  print("Normalizing lncRNA data...")
  lncrna_norm_matrix <- run_vst(lncrna_counts_filtered)
  
  print("--- mRNA and lncRNA normalization complete ---")
  print("Normalized mRNA matrix dimensions:")
  print(dim(mrna_norm_matrix))
  
  
  
  print("--- Checking for non-binary values in mo.data$mutation ---")
  print(table(mo.data$mutation))
  
  # A more concise way to see *only* the unique values:
  print("Unique values found:")
  print(unique(as.vector(mo.data$mutation)))
  
  print("--- Checking for NA values in mutation matrix ---")
  print(any(is.na(mo.data$mutation)))
  
  
  print("Fixing NA values in mutation matrix (setting NA to 0)...")
  
  # Get the matrix
  mut_matrix <- mo.data$mutation
  
  # Set any value that is NA to 0
  mut_matrix[is.na(mut_matrix)] <- 0
  
  # Put the fixed matrix back into the mo.data list
  mo.data$mutation <- mut_matrix
  
  # --- Verification ---
  print("--- Checking fixed matrix ---")
  print(paste("Any NAs left?", any(is.na(mo.data$mutation))))
  print("Unique values should now only be 0 and 1:")
  print(table(mo.data$mutation))
  print("--- Building the final 'mo.data' list for MOVICS ---")
  
  mo.data <- list(
    # Use the NEW normalized matrices
    mRNA = mrna_norm_matrix,
    lncRNA = lncrna_norm_matrix,
    
    # Use the ORIGINAL matrices, but filtered to the common patients
    mutation = binary_mutation_matrix[, common_patients],
    methylation = methylation_matrix_clean[, common_patients],
    cna = cna_matrix[, common_patients]
  )
  
  # --- Verification Step ---
  print("Final 'mo.data' list created. Checking dimensions (all should have 354 columns):")
  print(lapply(mo.data, function(x) dim(x)))
  
  
  
  
  # --- 1. Install/Load MOVICS library ---
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("MOVICS", quietly = TRUE)) {
    devtools::install_github("xlucpu/MOVICS")
  }
  library(MOVICS)
  
  args(getElites)
  
  # --- 1. Load the MOVICS library ---
  library(MOVICS)
  
  print("--- Step 1: Running Feature Selection (getElites) ---")
  
  # --- Continuous Data (mRNA, lncRNA, meth, cna) ---
  print("Selecting top 1000 for mRNA...")
  feat_mrna <- MOVICS::getElites(
    dat = mo.data$mRNA,
    method = "mad",
    elite.num = 1000  # <--- This is the correct argument
  )
  print(feat_mrna)
  
  print("Selecting top 1000 for lncRNA...")
  feat_lncrna <- MOVICS::getElites(
    dat = mo.data$lncRNA,
    method = "mad",
    elite.num = 1000  # <--- This is the correct argument
  )
  print(feat_lncrna)
  print("Selecting top 1000 for methylation...")
  feat_meth <- MOVICS::getElites(
    dat = mo.data$methylation,
    method = "mad",
    elite.num = 1000  # <--- This is the correct argument
  )
  print(feat_meth)
  print("Selecting top 1000 for CNA...")
  feat_cna <- MOVICS::getElites(
    dat = mo.data$cna,
    method = "mad",
    elite.num = 1000  # <--- This is the correct argument
  )
  
  # --- Binary Data (mutation) ---
  print("Selecting mutation features with > 10% frequency...")
  print("Selecting mutation features with > 10% frequency...")
  feat_mut <- MOVICS::getElites(
    dat = mo.data$mutation,
    method = "freq",
     # <--- This is the correct argument
    elite.num=1000
  )
  print(feat_mut)
  # --- 2. Combine results into the final 'feat.data' list ---
  input.data <- list(
    mRNA        = feat.data$mRNA$elite.dat,
    lncRNA      = feat.data$lncRNA$elite.dat,
    methylation = feat.data$methylation$elite.dat,
    cna         = feat.data$cna$elite.dat
    #his one was manual and is already correct
  )
  
  print("Selecting Top 1000 most mutated genes...")
  
  # (Make sure 'mo.data$mutation' is the clean, numeric version)
  # Calculate mutation *counts* for each gene (row)
  gene_counts <- rowSums(mo.data$mutation, na.rm = TRUE)
  
  # Get the names of the top 1000 genes
  top_1000_mut_genes <- names(sort(gene_counts, decreasing = TRUE)[1:1000])
  
  # Filter the matrix and add it to the list
  input.data$mutation <- mo.data$mutation[top_1000_mut_genes, ]
  
  # --- 4. Verification ---
  print("--- Final features selected per omics type: ---")
  print(lapply(input.data, dim))
  
  # --- 2. Verification ---
  # NOW, this command will work correctly.
  print("--- Final features selected per omics type: ---")
  print(lapply(input.data, dim))
  
  # --- Check the inputs to getElites ---
  
  print("--- Checking dimensions of matrices INSIDE mo.data ---")
  print(paste("mRNA dimensions:", paste(dim(mo.data$mRNA), collapse = " x ")))
  print(paste("lncRNA dimensions:", paste(dim(mo.data$lncRNA), collapse = " x ")))
  print(paste("Methylation dimensions:", paste(dim(mo.data$methylation), collapse = " x ")))
  print(paste("CNA dimensions:", paste(dim(mo.data$cna), collapse = " x ")))
  print(paste("Mutation dimensions:", paste(dim(mo.data$mutation), collapse = " x ")))
  
  
  # --- 1. Check ALL 5 matrices for NA values ---
  print("--- Checking for NA values ---")
  print(paste("NAs in mRNA:", any(is.na(mo.data$mRNA))))
  print(paste("NAs in lncRNA:", any(is.na(mo.data$lncRNA))))
  print(paste("NAs in methylation:", any(is.na(mo.data$methylation))))
  print(paste("NAs in CNA:", any(is.na(mo.data$cna))))
  print(paste("NAs in mutation:", any(is.na(mo.data$mutation))))
  
  # --- 2. Check the data type of the mutation matrix ---
  print("--- Checking mutation matrix data type ---")
  print(paste("Class of mutation matrix:", class(mo.data$mutation)))
  print(paste("Mode of mutation matrix:", mode(mo.data$mutation)))
  
  print("Fixing NA values in CNA matrix using k-NN imputation...")
  mo.data$cna <- impute.knn(mo.data$cna)$data 
  print("CNA matrix is clean.")
  
  # --- 1. Re-calculate the gene frequencies ---
  # (Make sure your 'mo.data$mutation' is the clean, numeric version)
  gene_frequencies <- rowSums(mo.data$mutation, na.rm = TRUE) / ncol(mo.data$mutation)
  
  # --- 2. Check the TOP 10 most mutated genes ---
  print("--- Top 10 most mutated genes (by frequency) ---")
  print(head(sort(gene_frequencies, decreasing = TRUE), 10))
  
  # --- 3. Check how many genes passed your filter ---
  genes_to_keep <- names(gene_frequencies[gene_frequencies > 0.10])
  print("--- Filter Check ---")
  print(paste("Total genes with > 10% mutation frequency:", length(genes_to_keep)))
  
  
  # --- 2. Determine Optimal Cluster Number (getClustNum) ---
  binary_vector <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  
  # --- 2. Run getClustNum ---
  print("--- Step 2: Calculating Optimal Cluster Number (getClustNum) ---")
  
  clust.num <- MOVICS::getClustNum(
    data = input.data,        # Your final data list
    is.binary = binary_vector,  # Tell the function which matrix is binary
    try.N.clust = 2:6           # Test k=2 through k=6 (from your plan)
    # We will use the defaults for center=TRUE and scale=TRUE
  )
  
  print("Optimal 'k' calculation complete. Check the Plots pane.")
  
  
  
  
  
  
  # --- 1. Define your 10 chosen algorithms ---
  my_algorithms <- c("CIMLR", "ConsensusClustering", "SNF", "iClusterBayes", 
                     "PINSPlus", "moCluster", "NEMO", "IntNMF", "COCA", "LRA")
  args(getMOIC)
  # --- 2. Define the *data type* for each matrix ---
  # (gaussian = continuous, binary = 0/1)
  # Our list order is: [mRNA, lncRNA, methylation, cna, mutation]
  data_type_vector <- c("gaussian", "gaussian", "gaussian", "gaussian", "binary")
  
  # --- 3. Run getMOIC with the CORRECT arguments ---
  print("--- Step 3: Running All 10 Clustering Algorithms (getMOIC) ---")
  print("This is the longest step and may take several minutes...")
  
  print("--- Step 3: Running All 10 Clustering Algorithms (getMOIC) ---")
  print("This is the longest step and may take several minutes...")
  
  # Define the data type for each matrix in 'input.data'
  # [mRNA, lncRNA, methylation, cna, mutation]
  my_algorithms <- c("CIMLR", "ConsensusClustering", "SNF", "iClusterBayes", 
                     "PINSPlus", "MoCluster", "NEMO", "IntNMF", "COCA", "LRAcluster")
  print("--- Step 3: Running All 10 Clustering Algorithms (getMOIC) ---")
  print("This is the longest step and may take several minutes...")
  data_type_vector <- c("gaussian", "gaussian", "gaussian", "gaussian", "binomial")
  # --- 3. Run getMOIC (with all arguments named) ---
  print("--- Step 3: Running All 10 Clustering Algorithms (getMOIC) ---")
  
  moic.res <- getMOIC(
    data = input.data,
    N.clust=2,
    methods = my_algorithms,
    type = data_type_vector,
 
  )
  
  print("All 10 clustering methods are complete.")
  
  # --- Step 5: Plot the Consensus Heatmap ---
  print("--- Step 5: Plotting Consensus Heatmap ---")
  
  plotConsensus(
    moic.res = moic.res,
    consensus.res = consensus.res,
    annot.col = c("MOIC_Subtype"), # Add the final C1/C2 labels
    seed = 123
  )
  
  
  # --- Step 8: Create the final clusters_df ---
  # This data frame is used for all downstream plots
  clusters_df <- data.frame(
    Patient = names(consensus.res$clust.res),
    MOIC_Subtype = paste0("C", consensus.res$clust.res)
  )
  rownames(clusters_df) <- clusters_df$Patient
  
  
  # --- A. MAFTOOLS ONCOPLOT ---
  # (This assumes 'lihc_maf_object_final' is in your environment)
  maf_clin_data <- lihc_maf_object_final@clinical.data
  maf_clin_data <- merge(maf_clin_data, clusters_df, by.x = "Tumor_Sample_Barcode", by.y = "Patient")
  rownames(maf_clin_data) <- maf_clin_data$Tumor_Sample_Barcode
  
  maf_subtypes <- read.maf(
    maf = lihc_maf_object_final@data,
    clinicalData = maf_clin_data,
    verbose = FALSE
  )
  
  oncoplot(
    maf = maf_subtypes,
    clinicalFeatures = "MOIC_Subtype", # This plots C1 vs C2
    sortByAnnotation = TRUE,
    top = 20
  )
  
  
  # --- B. LIMMA & GSEA PLOTS ---
  # (This assumes 'lihc_rna_se_final' is in your environment)
  rna_for_limma <- assay(lihc_rna_se_final)
  # ... (run limma as in your full script to get 'deg_results') ...
  # ... (run clusterProfiler as in your full script to get 'gse_go_results') ...
  #
  # dotplot(gse_go_results, showCategory = 20)
  
  
  # --- C. GSVA & IMMUNE BOXPLOTS ---
  # (This assumes 'mrna_matrix_symbol' is in your environment)
  # ... (run gsva/ssgsea as in your full script to get 'ssgsea_melt') ...
  #
  # ggplot(ssgsea_melt, aes(x = CellType, y = EnrichmentScore, fill = MOIC_Subtype)) +
  #   geom_boxplot() +
  #   stat_compare_means(aes(group = MOIC_Subtype), label = "p.signif"
  
  
  
  
  
  
  
  
  
  
  
  cluster1_raw <- c(
    "lihc_tcga:TCGA-2Y-A9GS-01", "lihc_tcga:TCGA-2Y-A9GT-01", "lihc_tcga:TCGA-2Y-A9GU-01",
    "lihc_tcga:TCGA-2Y-A9GV-01", "lihc_tcga:TCGA-2Y-A9GW-01", "lihc_tcga:TCGA-2Y-A9GX-01",
    "lihc_tcga:TCGA-2Y-A9GZ-01", "lihc_tcga:TCGA-2Y-A9H0-01", "lihc_tcga:TCGA-2Y-A9H4-01",
    "lihc_tcga:TCGA-2Y-A9H6-01", "lihc_tcga:TCGA-2Y-A9H9-01", "lihc_tcga:TCGA-2Y-A9HB-01",
    "lihc_tcga:TCGA-4R-AA8I-01", "lihc_tcga:TCGA-5C-AAPD-01", "lihc_tcga:TCGA-5R-AA1D-01",
    "lihc_tcga:TCGA-5R-AAAM-01", "lihc_tcga:TCGA-BC-A10S-01", "lihc_tcga:TCGA-BC-A10T-01",
    "lihc_tcga:TCGA-BC-A10X-01", "lihc_tcga:TCGA-BC-A10Y-01", "lihc_tcga:TCGA-BC-A10Z-01",
    "lihc_tcga:TCGA-BC-A110-01", "lihc_tcga:TCGA-BC-A216-01", "lihc_tcga:TCGA-BC-A217-01",
    "lihc_tcga:TCGA-BC-A3KF-01", "lihc_tcga:TCGA-BC-A3KG-01", "lihc_tcga:TCGA-BC-A5W4-01",
    "lihc_tcga:TCGA-BC-A8YO-01", "lihc_tcga:TCGA-BD-A2L6-01", "lihc_tcga:TCGA-BD-A3EP-01",
    "lihc_tcga:TCGA-BD-A3ER-01", "lihc_tcga:TCGA-BW-A5NO-01", "lihc_tcga:TCGA-CC-5259-01",
    "lihc_tcga:TCGA-CC-5262-01", "lihc_tcga:TCGA-CC-5264-01", "lihc_tcga:TCGA-CC-A3MB-01",
    "lihc_tcga:TCGA-CC-A5UE-01", "lihc_tcga:TCGA-CC-A7IF-01", "lihc_tcga:TCGA-CC-A8HT-01",
    "lihc_tcga:TCGA-CC-A8HU-01", "lihc_tcga:TCGA-CC-A9FS-01", "lihc_tcga:TCGA-DD-A115-01",
    "lihc_tcga:TCGA-DD-A119-01", "lihc_tcga:TCGA-DD-A11A-01", "lihc_tcga:TCGA-DD-A11B-01",
    "lihc_tcga:TCGA-DD-A11C-01", "lihc_tcga:TCGA-DD-A11D-01", "lihc_tcga:TCGA-DD-A1EB-01",
    "lihc_tcga:TCGA-DD-A1EG-01", "lihc_tcga:TCGA-DD-A1EJ-01", "lihc_tcga:TCGA-DD-A1EK-01",
    "lihc_tcga:TCGA-DD-A39X-01", "lihc_tcga:TCGA-DD-A39Y-01", "lihc_tcga:TCGA-DD-A39Z-01",
    "lihc_tcga:TCGA-DD-A3A1-01", "lihc_tcga:TCGA-DD-A3A2-01", "lihc_tcga:TCGA-DD-A3A4-01",
    "lihc_tcga:TCGA-DD-A3A5-01", "lihc_tcga:TCGA-DD-A3A7-01", "lihc_tcga:TCGA-DD-A4NB-01",
    "lihc_tcga:TCGA-DD-A4NE-01", "lihc_tcga:TCGA-DD-A4NI-01", "lihc_tcga:TCGA-DD-A4NJ-01",
    "lihc_tcga:TCGA-DD-A4NK-01", "lihc_tcga:TCGA-DD-A4NL-01", "lihc_tcga:TCGA-DD-A4NO-01",
    "lihc_tcga:TCGA-DD-A4NQ-01", "lihc_tcga:TCGA-DD-A4NS-01", "lihc_tcga:TCGA-DD-A4NV-01",
    "lihc_tcga:TCGA-DD-A73A-01", "lihc_tcga:TCGA-DD-A73B-01", "lihc_tcga:TCGA-DD-A73C-01",
    "lihc_tcga:TCGA-DD-A73F-01", "lihc_tcga:TCGA-DD-AAC9-01", "lihc_tcga:TCGA-DD-AACE-01",
    "lihc_tcga:TCGA-DD-AACF-01", "lihc_tcga:TCGA-DD-AACG-01", "lihc_tcga:TCGA-DD-AACI-01",
    "lihc_tcga:TCGA-DD-AACO-01", "lihc_tcga:TCGA-DD-AACQ-01", "lihc_tcga:TCGA-DD-AACS-01",
    "lihc_tcga:TCGA-DD-AACT-01", "lihc_tcga:TCGA-DD-AACU-01", "lihc_tcga:TCGA-DD-AACV-01",
    "lihc_tcga:TCGA-DD-AACW-01", "lihc_tcga:TCGA-DD-AACY-01", "lihc_tcga:TCGA-DD-AAD2-01",
    "lihc_tcga:TCGA-DD-AAD3-01", "lihc_tcga:TCGA-DD-AAD8-01", "lihc_tcga:TCGA-DD-AADF-01",
    "lihc_tcga:TCGA-DD-AADJ-01", "lihc_tcga:TCGA-DD-AADL-01", "lihc_tcga:TCGA-DD-AADM-01",
    "lihc_tcga:TCGA-DD-AADN-01", "lihc_tcga:TCGA-DD-AADP-01", "lihc_tcga:TCGA-DD-AADR-01",
    "lihc_tcga:TCGA-DD-AADS-01", "lihc_tcga:TCGA-DD-AADV-01", "lihc_tcga:TCGA-DD-AAE1-01",
    "lihc_tcga:TCGA-DD-AAE2-01", "lihc_tcga:TCGA-DD-AAE3-01", "lihc_tcga:TCGA-DD-AAE4-01",
    "lihc_tcga:TCGA-DD-AAE6-01", "lihc_tcga:TCGA-DD-AAE7-01", "lihc_tcga:TCGA-DD-AAEA-01",
    "lihc_tcga:TCGA-DD-AAED-01", "lihc_tcga:TCGA-DD-AAEG-01", "lihc_tcga:TCGA-DD-AAEH-01",
    "lihc_tcga:TCGA-DD-AAEI-01", "lihc_tcga:TCGA-DD-AAVU-01", "lihc_tcga:TCGA-DD-AAVW-01",
    "lihc_tcga:TCGA-DD-AAVZ-01", "lihc_tcga:TCGA-DD-AAW0-01", "lihc_tcga:TCGA-ED-A459-01",
    "lihc_tcga:TCGA-ED-A4XI-01", "lihc_tcga:TCGA-ED-A627-01", "lihc_tcga:TCGA-ED-A7XO-01",
    "lihc_tcga:TCGA-EP-A2KB-01", "lihc_tcga:TCGA-EP-A3JL-01", "lihc_tcga:TCGA-ES-A2HS-01",
    "lihc_tcga:TCGA-ES-A2HT-01", "lihc_tcga:TCGA-FV-A23B-01", "lihc_tcga:TCGA-FV-A2QQ-01",
    "lihc_tcga:TCGA-FV-A2QR-01", "lihc_tcga:TCGA-FV-A3I1-01", "lihc_tcga:TCGA-FV-A3R2-01",
    "lihc_tcga:TCGA-FV-A495-01", "lihc_tcga:TCGA-FV-A496-01", "lihc_tcga:TCGA-G3-A25S-01",
    "lihc_tcga:TCGA-G3-A25U-01", "lihc_tcga:TCGA-G3-A25V-01", "lihc_tcga:TCGA-G3-A25Z-01",
    "lihc_tcga:TCGA-G3-A3CG-01", "lihc_tcga:TCGA-G3-A3CH-01", "lihc_tcga:TCGA-G3-A3CI-01",
    "lihc_tcga:TCGA-G3-A5SI-01", "lihc_tcga:TCGA-G3-A5SK-01", "lihc_tcga:TCGA-G3-A5SM-01",
    "lihc_tcga:TCGA-G3-A7M5-01", "lihc_tcga:TCGA-G3-A7M7-01", "lihc_tcga:TCGA-G3-A7M8-01",
    "lihc_tcga:TCGA-G3-AAUZ-01", "lihc_tcga:TCGA-G3-AAV0-01", "lihc_tcga:TCGA-GJ-A6C0-01",
    "lihc_tcga:TCGA-GJ-A9DB-01", "lihc_tcga:TCGA-HP-A5MZ-01", "lihc_tcga:TCGA-HP-A5N0-01",
    "lihc_tcga:TCGA-K7-A5RF-01", "lihc_tcga:TCGA-K7-A5RG-01", "lihc_tcga:TCGA-K7-A6G5-01",
    "lihc_tcga:TCGA-KR-A7K2-01", "lihc_tcga:TCGA-KR-A7K8-01", "lihc_tcga:TCGA-LG-A6GG-01",
    "lihc_tcga:TCGA-LG-A9QD-01", "lihc_tcga:TCGA-MI-A75E-01", "lihc_tcga:TCGA-MI-A75G-01",
    "lihc_tcga:TCGA-MI-A75I-01", "lihc_tcga:TCGA-MR-A520-01", "lihc_tcga:TCGA-NI-A4U2-01",
    "lihc_tcga:TCGA-O8-A75V-01", "lihc_tcga:TCGA-PD-A5DF-01", "lihc_tcga:TCGA-QA-A7B7-01",
    "lihc_tcga:TCGA-RC-A6M5-01", "lihc_tcga:TCGA-RC-A7SB-01", "lihc_tcga:TCGA-RC-A7SF-01",
    "lihc_tcga:TCGA-RC-A7SK-01", "lihc_tcga:TCGA-T1-A6J8-01", "lihc_tcga:TCGA-UB-A7MB-01",
    "lihc_tcga:TCGA-UB-AA0U-01", "lihc_tcga:TCGA-UB-AA0V-01", "lihc_tcga:TCGA-WJ-A86L-01",
    "lihc_tcga:TCGA-WQ-AB4B-01", "lihc_tcga:TCGA-WX-AA46-01", "lihc_tcga:TCGA-XR-A8TC-01",
    "lihc_tcga:TCGA-XR-A8TD-01", "lihc_tcga:TCGA-XR-A8TG-01", "lihc_tcga:TCGA-ZP-A9CV-01",
    "lihc_tcga:TCGA-ZP-A9CY-01", "lihc_tcga:TCGA-ZP-A9CZ-01", "lihc_tcga:TCGA-ZP-A9D0-01",
    "lihc_tcga:TCGA-ZP-A9D1-01", "lihc_tcga:TCGA-ZS-A9CD-01", "lihc_tcga:TCGA-ZS-A9CE-01",
    "lihc_tcga:TCGA-ZS-A9CF-01", "lihc_tcga:TCGA-ZS-A9CF-02"
  )
  
  cluster2_raw <- c(
    "lihc_tcga:TCGA-2V-A95S-01", "lihc_tcga:TCGA-2Y-A9GY-01", "lihc_tcga:TCGA-2Y-A9H2-01",
    "lihc_tcga:TCGA-2Y-A9H3-01", "lihc_tcga:TCGA-2Y-A9H5-01", "lihc_tcga:TCGA-2Y-A9H7-01",
    "lihc_tcga:TCGA-2Y-A9H8-01", "lihc_tcga:TCGA-5C-A9VG-01", "lihc_tcga:TCGA-5C-A9VH-01",
    "lihc_tcga:TCGA-BC-4072-01", "lihc_tcga:TCGA-BC-4073-01", "lihc_tcga:TCGA-BC-A10Q-01",
    "lihc_tcga:TCGA-BC-A10R-01", "lihc_tcga:TCGA-BC-A10W-01", "lihc_tcga:TCGA-BC-A112-01",
    "lihc_tcga:TCGA-BC-A69H-01", "lihc_tcga:TCGA-BW-A5NP-01", "lihc_tcga:TCGA-CC-5258-01",
    "lihc_tcga:TCGA-CC-5260-01", "lihc_tcga:TCGA-CC-5261-01", "lihc_tcga:TCGA-CC-5263-01",
    "lihc_tcga:TCGA-CC-A123-01", "lihc_tcga:TCGA-CC-A1HT-01", "lihc_tcga:TCGA-CC-A3M9-01",
    "lihc_tcga:TCGA-CC-A3MA-01", "lihc_tcga:TCGA-CC-A3MC-01", "lihc_tcga:TCGA-CC-A5UC-01",
    "lihc_tcga:TCGA-CC-A5UD-01", "lihc_tcga:TCGA-CC-A7IE-01", "lihc_tcga:TCGA-CC-A7IG-01",
    "lihc_tcga:TCGA-CC-A7II-01", "lihc_tcga:TCGA-CC-A7IJ-01", "lihc_tcga:TCGA-CC-A8HS-01",
    "lihc_tcga:TCGA-CC-A8HV-01", "lihc_tcga:TCGA-CC-A9FU-01", "lihc_tcga:TCGA-CC-A9FV-01",
    "lihc_tcga:TCGA-DD-A113-01", "lihc_tcga:TCGA-DD-A114-01", "lihc_tcga:TCGA-DD-A118-01",
    "lihc_tcga:TCGA-DD-A1EC-01", "lihc_tcga:TCGA-DD-A1EF-01", "lihc_tcga:TCGA-DD-A1EH-01",
    "lihc_tcga:TCGA-DD-A1EI-01", "lihc_tcga:TCGA-DD-A3A6-01", "lihc_tcga:TCGA-DD-A3A9-01",
    "lihc_tcga:TCGA-DD-A4NA-01", "lihc_tcga:TCGA-DD-A4ND-01", "lihc_tcga:TCGA-DD-A4NH-01",
    "lihc_tcga:TCGA-DD-A4NN-01", "lihc_tcga:TCGA-DD-A4NR-01", "lihc_tcga:TCGA-DD-A73G-01",
    "lihc_tcga:TCGA-DD-AA3A-01", "lihc_tcga:TCGA-DD-AACC-01", "lihc_tcga:TCGA-DD-AACH-01",
    "lihc_tcga:TCGA-DD-AACL-01", "lihc_tcga:TCGA-DD-AACN-01", "lihc_tcga:TCGA-DD-AACP-01",
    "lihc_tcga:TCGA-DD-AACZ-01", "lihc_tcga:TCGA-DD-AAD1-01", "lihc_tcga:TCGA-DD-AAD5-01",
    "lihc_tcga:TCGA-DD-AADA-01", "lihc_tcga:TCGA-DD-AADB-01", "lihc_tcga:TCGA-DD-AADC-01",
    "lihc_tcga:TCGA-DD-AADI-01", "lihc_tcga:TCGA-DD-AADK-01", "lihc_tcga:TCGA-DD-AADO-01",
    "lihc_tcga:TCGA-DD-AADW-01", "lihc_tcga:TCGA-DD-AADY-01", "lihc_tcga:TCGA-DD-AAE0-01",
    "lihc_tcga:TCGA-DD-AAVQ-01", "lihc_tcga:TCGA-DD-AAVR-01", "lihc_tcga:TCGA-DD-AAVS-01",
    "lihc_tcga:TCGA-DD-AAVV-01", "lihc_tcga:TCGA-ED-A5KG-01", "lihc_tcga:TCGA-ED-A66X-01",
    "lihc_tcga:TCGA-ED-A66Y-01", "lihc_tcga:TCGA-ED-A7PX-01", "lihc_tcga:TCGA-ED-A7PY-01",
    "lihc_tcga:TCGA-ED-A82E-01", "lihc_tcga:TCGA-ED-A8O5-01", "lihc_tcga:TCGA-ED-A8O6-01",
    "lihc_tcga:TCGA-ED-A97K-01", "lihc_tcga:TCGA-EP-A2KA-01", "lihc_tcga:TCGA-FV-A3I0-01",
    "lihc_tcga:TCGA-FV-A3R3-01", "lihc_tcga:TCGA-FV-A4ZP-01", "lihc_tcga:TCGA-G3-A25T-01",
    "lihc_tcga:TCGA-G3-A25X-01", "lihc_tcga:TCGA-G3-A25Y-01", "lihc_tcga:TCGA-G3-A5SJ-01",
    "lihc_tcga:TCGA-G3-A7M6-01", "lihc_tcga:TCGA-G3-A7M9-01", "lihc_tcga:TCGA-G3-AAV6-01",
    "lihc_tcga:TCGA-G3-AAV7-01", "lihc_tcga:TCGA-GJ-A3OU-01", "lihc_tcga:TCGA-K7-AAU7-01",
    "lihc_tcga:TCGA-MR-A8JO-01", "lihc_tcga:TCGA-RC-A6M3-01", "lihc_tcga:TCGA-RC-A6M6-01",
    "lihc_tcga:TCGA-RC-A7S9-01", "lihc_tcga:TCGA-RC-A7SH-01", "lihc_tcga:TCGA-UB-A7MA-01",
    "lihc_tcga:TCGA-UB-A7ME-01", "lihc_tcga:TCGA-UB-A7MF-01", "lihc_tcga:TCGA-WQ-A9G7-01",
    "lihc_tcga:TCGA-WX-AA44-01", "lihc_tcga:TCGA-YA-A8S7-01", "lihc_tcga:TCGA-ZP-A9D2-01"
  )
  
  cluster3_raw <- c(
    "lihc_tcga:TCGA-2Y-A9H1-01", "lihc_tcga:TCGA-2Y-A9HA-01", "lihc_tcga:TCGA-3K-AAZ8-01",
    "lihc_tcga:TCGA-5R-AA1C-01", "lihc_tcga:TCGA-BC-A10U-01", "lihc_tcga:TCGA-BC-A69I-01",
    "lihc_tcga:TCGA-BW-A5NQ-01", "lihc_tcga:TCGA-CC-A7IH-01", "lihc_tcga:TCGA-CC-A7IK-01",
    "lihc_tcga:TCGA-CC-A7IL-01", "lihc_tcga:TCGA-CC-A9FW-01", "lihc_tcga:TCGA-DD-A116-01",
    "lihc_tcga:TCGA-DD-A1EA-01", "lihc_tcga:TCGA-DD-A1ED-01", "lihc_tcga:TCGA-DD-A1EE-01",
    "lihc_tcga:TCGA-DD-A1EL-01", "lihc_tcga:TCGA-DD-A39V-01", "lihc_tcga:TCGA-DD-A39W-01",
    "lihc_tcga:TCGA-DD-A3A3-01", "lihc_tcga:TCGA-DD-A3A8-01", "lihc_tcga:TCGA-DD-A4NF-01",
    "lihc_tcga:TCGA-DD-A4NG-01", "lihc_tcga:TCGA-DD-A4NP-01", "lihc_tcga:TCGA-DD-A73D-01",
    "lihc_tcga:TCGA-DD-A73E-01", "lihc_tcga:TCGA-DD-AAC8-01", "lihc_tcga:TCGA-DD-AACA-01",
    "lihc_tcga:TCGA-DD-AACA-02", "lihc_tcga:TCGA-DD-AACB-01", "lihc_tcga:TCGA-DD-AACD-01",
    "lihc_tcga:TCGA-DD-AACJ-01", "lihc_tcga:TCGA-DD-AACK-01", "lihc_tcga:TCGA-DD-AACX-01",
    "lihc_tcga:TCGA-DD-AAD0-01", "lihc_tcga:TCGA-DD-AAD6-01", "lihc_tcga:TCGA-DD-AADD-01",
    "lihc_tcga:TCGA-DD-AADG-01", "lihc_tcga:TCGA-DD-AADQ-01", "lihc_tcga:TCGA-DD-AADU-01",
    "lihc_tcga:TCGA-DD-AAE9-01", "lihc_tcga:TCGA-DD-AAEB-01", "lihc_tcga:TCGA-DD-AAEE-01",
    "lihc_tcga:TCGA-DD-AAEK-01", "lihc_tcga:TCGA-DD-AAVP-01", "lihc_tcga:TCGA-DD-AAVX-01",
    "lihc_tcga:TCGA-DD-AAVY-01", "lihc_tcga:TCGA-DD-AAW1-01", "lihc_tcga:TCGA-DD-AAW2-01",
    "lihc_tcga:TCGA-DD-AAW3-01", "lihc_tcga:TCGA-ED-A7PZ-01", "lihc_tcga:TCGA-ED-A7XP-01",
    "lihc_tcga:TCGA-EP-A12J-01", "lihc_tcga:TCGA-EP-A26S-01", "lihc_tcga:TCGA-EP-A2KC-01",
    "lihc_tcga:TCGA-EP-A3RK-01", "lihc_tcga:TCGA-FV-A4ZQ-01", "lihc_tcga:TCGA-G3-A3CJ-01",
    "lihc_tcga:TCGA-G3-A3CK-01", "lihc_tcga:TCGA-G3-A5SL-01", "lihc_tcga:TCGA-G3-A6UC-01",
    "lihc_tcga:TCGA-G3-AAV1-01", "lihc_tcga:TCGA-G3-AAV2-01", "lihc_tcga:TCGA-G3-AAV3-01",
    "lihc_tcga:TCGA-G3-AAV4-01", "lihc_tcga:TCGA-G3-AAV5-01", "lihc_tcga:TCGA-KR-A7K0-01",
    "lihc_tcga:TCGA-KR-A7K7-01", "lihc_tcga:TCGA-LG-A9QC-01", "lihc_tcga:TCGA-MI-A75C-01",
    "lihc_tcga:TCGA-MI-A75H-01", "lihc_tcga:TCGA-NI-A8LF-01", "lihc_tcga:TCGA-RC-A6M4-01",
    "lihc_tcga:TCGA-RG-A7D4-01", "lihc_tcga:TCGA-UB-A7MC-01", "lihc_tcga:TCGA-UB-A7MD-01",
    "lihc_tcga:TCGA-WX-AA47-01", "lihc_tcga:TCGA-XR-A8TE-01", "lihc_tcga:TCGA-XR-A8TF-01",
    "lihc_tcga:TCGA-ZP-A9D4-01", "lihc_tcga:TCGA-ZS-A9CG-01."
  )
  
  # Check the counts to verify data entry
  # Expected: 185, 108, 80
  print("--- Original Counts ---")
  print(paste("Cluster 1:", length(cluster1_raw)))
  print(paste("Cluster 2:", length(cluster2_raw)))
  print(paste("Cluster 3:", length(cluster3_raw)))
  
  
  # --- 2. Extraction Methods ---
  
  # We can put all the raw data into a list to process them all at once
  all_clusters_raw <- list(
    cluster1 = cluster1_raw,
    cluster2 = cluster2_raw,
    cluster3 = cluster3_raw
  )
  
  # --- Method 1: Base R using sub() and substr() ---
  # This is a simple, two-step approach that is easy to understand.
  # 1. sub(".*:", "", x) -> Removes everything before and including the colon
  # 2. substr(..., 1, 12) -> Takes the first 12 characters of the result
  
  extract_ids_simple <- function(id_vector) {
    # Remove the prefix (e.g., "lihc_tcga:")
    ids_no_prefix <- sub(".*:", "", id_vector)
    # Extract the first 12 characters
    patient_ids <- substr(ids_no_prefix, 1, 12)
    return(patient_ids)
  }
  
  # Apply this function to all clusters
  all_clusters_simple <- lapply(all_clusters_raw, extract_ids_simple)
  
  cat("\n--- Results from Method 1 (sub + substr) ---\n")
  print(lapply(all_clusters_simple, head))
  
  
  # --- Method 2: Base R using regexpr() and regmatches() ---
  # This is more precise, as it looks for the specific TCGA patient ID format.
  # Pattern: "TCGA-" + 2 alphanumeric chars + "-" + 4 alphanumeric chars
  # This is generally the recommended base R method.
  
  extract_ids_regex <- function(id_vector) {
    # Define the regex pattern for the 12-character patient ID
    pattern <- "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
    
    # Find all matches
    matches <- regexpr(pattern, id_vector)
    
    # Extract the matched strings
    patient_ids <- regmatches(id_vector, matches)
    
    return(patient_ids)
  }
  
  # Apply this function to all clusters
  all_clusters_regex <- lapply(all_clusters_raw, extract_ids_regex)
  
  cat("\n--- Results from Method 2 (Base R Regex) ---\n")
  print(lapply(all_clusters_regex, head))
  
  # Check counts
  cat("\n--- Counts from Method 2 ---\n")
  print(lapply(all_clusters_regex, length))
  
  
  # --- Method 3: Using the 'stringr' package (Tidyverse) ---
  # This is often the cleanest way if you use the Tidyverse.
  # You may need to install it first: install.packages("stringr")
  
  # We create a function, but you could also just use the stringr::str_extract line
  extract_ids_stringr <- function(id_vector) {
    # Check if stringr is loaded, if not, try to load it
    if (!requireNamespace("stringr", quietly = TRUE)) {
      message("stringr package not found. Please install it using: install.packages(\"stringr\")")
      return(NULL)
    }
    
    # Define the same pattern
    pattern <- "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
    
    # Extract the first match
    patient_ids <- stringr::str_extract(id_vector, pattern)
    
    # str_extract returns NA for no match, which is fine
    return(patient_ids)
  }
  
  # Apply this function to all clusters
  all_clusters_stringr <- lapply(all_clusters_raw, extract_ids_stringr)
  
  cat("\n--- Results from Method 3 (stringr) ---\n")
  # We add na.rm = TRUE to head() in case the package isn't installed
  print(lapply(all_clusters_stringr, head, na.rm = TRUE))
  
  
  # --- 4. Final Verification ---
  # Let's check the special cases, like the -02 samples and the one with a period
  cat("\n--- Special Case Checks (using Method 2) ---\n")
  
  # Check the sample with a trailing period
  cat("Original (Cluster 3 tail):", tail(cluster3_raw, 1), "\n")
  cat("Extracted (Cluster 3 tail):", tail(all_clusters_regex$cluster3, 1), "\n\n")
  
  # Check the -02 sample
  cat("Original (Cluster 1 tail):", tail(cluster1_raw, 1), "\n")
  cat("Extracted (Cluster 1 tail):", tail(all_clusters_regex$cluster1, 1), "\n")
  
  
  # You can access the final, clean vectors of IDs like this:
  # (Assuming you chose Method 2)
  cluster1_ids <- all_clusters_regex$cluster1
  cluster2_ids <- all_clusters_regex$cluster2
  cluster3_ids <- all_clusters_regex$cluster3
  
  cat("\n--- Clean Cluster 1 IDs ---\n")
  cat(cluster1_ids, sep = "\n")
  
  cat("\n--- Clean Cluster 2 IDs ---\n")
  cat(cluster2_ids, sep = "\n")
  
  cat("\n--- Clean Cluster 3 IDs ---\n")
  cat(cluster3_ids, sep = "\n")
  

  all_patient_ids <- c(cluster1_ids, cluster2_ids, cluster3_ids)
  
  # Print the total count to verify
  # Expected: 185 + 108 + 80 = 373
  cat("Total number of patient IDs (all clusters):", length(all_patient_ids), "\n")
  
  # Print the first few IDs from the combined list
  cat("\nHead of combined list:\n")
  print(head(all_patient_ids))
  
  # Print the last few IDs from the combined list
  cat("\nTail of combined list:\n")
  print(tail(all_patient_ids))
  
  
  
  # ------------------------------------------------------------
  # 1. Install & load library
  # ------------------------------------------------------------
  if (!require(UCSCXenaTools)) install.packages("UCSCXenaTools")
  library(UCSCXenaTools)
  
  # ------------------------------------------------------------
  # 2. Specify dataset (TCGA-LIHC miRNA stem-loop expression)
  # ------------------------------------------------------------
  dataset <- "TCGA-LIHC.mirna.tsv"
  hub <- "https://gdc.xenahubs.net"
  
  # ------------------------------------------------------------
  # 3. Query + Download metadata and data
  # ------------------------------------------------------------
  xe = XenaGenerate() %>%
    XenaFilter(filterDatasets = dataset) %>%
    XenaFilter(filterHosts = hub) %>%
    XenaQuery() %>%
    XenaDownload()
  
  # ------------------------------------------------------------
  # 4. Load matrix into R
  # ------------------------------------------------------------
  library(data.table)
  
  # Locate the downloaded file (Xena always places it under temp)
  filepath <- list.files(
    path = tempdir(),
    pattern = "TCGA-LIHC\\.mirna\\.tsv\\.gz$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Read the matrix
  miRNA_df <- fread(filepath, sep = "\t", header = TRUE)
  # First column is miRNA IDs → keep it unchanged
  colnames(miRNA_df)[-1] <- substr(colnames(miRNA_df)[-1], 1, 12)
  
  # Inspect
  head(colnames(miRNA_df))
  
  # Inspect
  dim(miRNA_df)
  head(miRNA_df[, 1:5])
  
  miRNA_mat <- as.data.frame(miRNA_df)
  rownames(miRNA_mat) <- miRNA_mat[[1]]
  miRNA_mat[[1]] <- NULL
  head(miRNA_mat[, 1:5])
  # Convert to matrix
  miRNA_mat <- as.matrix(miRNA_mat)
  
  patients_to_keep <- intersect(colnames(miRNA_mat), all_patient_ids)
  
  # Print a summary of the overlap
  cat(paste("Found", length(patients_to_keep), "of your", length(all_patient_ids), "patients in the miRNA data.\n"))
  
  # Subset the matrix to *only* these patients (columns)
  miRNA_mat_filtered <- miRNA_mat[, patients_to_keep]
  
  # Inspect new dimensions
  # Should have the same number of rows (miRNAs) but fewer columns (patients)
  dim(miRNA_mat)
  dim(miRNA_mat_filtered)
  
  
  # ------------------------------------------------------------
  # 4. Select top variable miRNAs (top 100)
  # ------------------------------------------------------------
  
  # Compute variance across samples (rows = miRNAs)
  vars <- apply(miRNA_mat_filtered, 1, var)
  
  # Choose top 100 most variable miRNAs
  top_miRNA <- names(sort(vars, decreasing = TRUE))[1:500]
  
  # Subset matrix
  mat_top <- miRNA_mat_filtered[top_miRNA, ]
  
  # Get all the patient IDs that are actually in your final matrix
  all_patients_in_matrix <- colnames(mat_top)
  
  # 1. Create a data frame with patient IDs as row names
  annotation_df <- data.frame(
    Cluster = rep(NA, length(all_patients_in_matrix)),
    row.names = all_patients_in_matrix
  )
  
  # 2. Find which of your cluster patients are in the matrix
  # (We use intersect to be safe)
  c1_in_mat <- intersect(all_patients_in_matrix, cluster1_ids)
  c2_in_mat <- intersect(all_patients_in_matrix, cluster2_ids)
  c3_in_mat <- intersect(all_patients_in_matrix, cluster3_ids)
  
  # 3. Fill in the 'Cluster' column based on your lists
  annotation_df[c1_in_mat, "Cluster"] <- "Cluster 1"
  annotation_df[c2_in_mat, "Cluster"] <- "Cluster 2"
  annotation_df[c3_in_mat, "Cluster"] <- "Cluster 3"
  
  # Check the first few rows
  head(annotation_df)
  print(table(annotation_df$Cluster))
  
  
  # --- 1. Normalize data: 0-1 range for each row (miRNA) ---
  
 
  
  # Apply the function to each row (MARGIN = 1)
  # We must transpose the result (t()) because apply returns results in columns
  mat_top_norm_01 <- t(scale(t(mat_top)))
  
  # Verify the result (all row mins should be 0, all row maxs should be 1)
  cat("\n--- Verifying 0-1 Normalization ---\n")
  print(paste("Row Mins (all should be 0):", head(apply(mat_top_norm_01, 1, min))))
  print(paste("Row Maxs (all should be 1):", head(apply(mat_top_norm_01, 1, max))))
  
  
  # This is the 500 miRNA (rows) x 373 patient (cols) matrix, 0-1 scaled
  # We transpose it to be 373 patients (rows) x 500 miRNA (cols)
  data_for_kmeans <- t(mat_top_norm_01) 
  
  set.seed(123) # For reproducible results
  
  # --- Run for k=2 ---
  kmeans_k2 <- kmeans(data_for_kmeans, centers = 2, nstart = 25)
  
  # --- Run for k=3 ---
  kmeans_k3 <- kmeans(data_for_kmeans, centers = 3, nstart = 25)
  
  # This assumes 'annotation_df' already exists with your original 'Cluster' column
  
  # Add k=2 results as a new column
  annotation_df$KMeans_k2 <- as.factor(kmeans_k2$cluster)
  
  # Add k=3 results as a new column
  # (This replaces the old 'KMeans_Cluster' if it existed, which is fine)
  annotation_df$KMeans_k3 <- as.factor(kmeans_k3$cluster)
  
  # Check the updated annotation (it now has 3+ columns)
  cat("\n--- Annotation updated with k=2 and k=3 ---\n")
  print(head(annotation_df))
  # ------------------------------------------------------------
  # 5. Generate Heatmap
  # ------------------------------------------------------------
  
  # Option 1: Unsupervised clustering (see if your groups match)
  # You may need this package for the color palette
  # install.packages("RColorBrewer")
  
  # 1. Define a sequential color palette
  # This creates 100 colors from Yellow to Orange to Red
  my_colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))(100)
  
  # 1. Get the order of patients based on k=3 results
  col_order_nested <- order(
   
    annotation_df$Cluster,
    
    annotation_df$KMeans_k3,
    annotation_df$KMeans_k2
  )
  # 2. Re-order the matrix and annotation
  mat_top_ordered <- mat_top_norm_01[, col_order_nested]
  annotation_df_ordered <- annotation_df[col_order_nested, , drop = FALSE]
  # 3. Plot, turning OFF column clustering
  pheatmap(
    mat_top_ordered,         # <-- Use the re-ordered matrix
    scale = "none",
    color = my_colors,
    annotation_col = annotation_df_ordered, # <-- Use the re-ordered annotation
    cluster_cols = FALSE,          # <-- Turn off H-clustering for columns
    cluster_rows = TRUE,           # <-- Keep clustering miRNAs
    fontsize_row = 6,
    fontsize_col = 6,
    show_colnames = FALSE,
    main = "TCGA-LIHC: Grouped by K-Means (k=3)"
  )
  
  
  library(factoextra)
  
  # --- 2. Define the Data for Scoring ---
  # This is the 500 miRNA (rows) x ~373 patient (cols) matrix, 0-1 scaled
  # We transpose it to be ~373 patients (rows) x 500 miRNA (cols)
  data_for_kmeans <- t(mat_top_norm_01)
  
  # Set seed for reproducible results
  set.seed(123) 
  
  # --- 3. Compute Scores: Elbow Method (WSS) ---
  # This method plots the total within-cluster sum of squares (WSS).
  # You look for a "bend" or "elbow" in the plot.
  # After the elbow, adding more clusters doesn't explain much more variance.
  
  cat("\n--- Generating Elbow Method Plot (k=1 to 6) ---\n")
  fviz_nbclust(
    data_for_kmeans, 
    kmeans, 
    method = "wss",  # "wss" stands for Within Sum of Squares
    k.max = 6,         # Test up to k=6
    nstart = 25        # Run kmeans 25 times per k to find a good solution
  )
  
  # --- 4. Compute Scores: Silhouette Method (Recommended) ---
  # This method measures how similar a patient is to its own cluster
  # compared to other clusters. A higher average silhouette width (closer to 1) is better.
  # **You look for the highest peak.**
  
  cat("\n--- Generating Silhouette Method Plot (k=2 to 6) ---\n")
  fviz_nbclust(
    data_for_kmeans, 
    kmeans, 
    method = "silhouette",
    k.max = 6,         # Test up to k=6
    nstart = 25
  )
  
  
  
  
  
  # --- 1. Install and Load Packages ---
  # We only need these two packages for this script
  if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap")
  
  if (!requireNamespace("clValid", quietly = TRUE))
    install.packages("clValid")
  if (!requireNamespace("factoextra", quietly = TRUE))
    install.packages("factoextra")
  
  # Load libraries
  library(pheatmap)
  library(clValid)

  
  # --- 0. Load Required Libraries ---
  # We need 'cluster' for pam() and 'dplyr' for arrange()
  library(cluster)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  
  # --- Helper function to find where to draw lines between clusters ---
  find_gaps <- function(annotation_df, variable) {
    # Sort the annotation frame by the variable
    annotation_sorted <- annotation_df[order(annotation_df[[variable]]), , drop = FALSE]
    # Find the gaps
    gaps <- which(annotation_sorted[[variable]][-1] != annotation_sorted[[variable]][-nrow(annotation_sorted)])
    return(gaps)
  }
  
  # --- 1. Load mRNA Data (TCGA-LIHC.star_counts.tsv.gz) ---
  # (Assuming 'all_patient_ids', 'cluster1_ids', etc. are in your environment)
  
  cat("Downloading TCGA-LIHC.star_counts.tsv.gz...\n")
  url_mrna <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz"
  dest_file_mrna <- "TCGA-LIHC.star_counts.tsv.gz"
  
  if (!file.exists(dest_file_mrna)) {
    download.file(url_mrna, dest_file_mrna)
  }
  
  mrna_df <- read.delim(gzfile(dest_file_mrna), check.names = FALSE)
  cat("mRNA Download and load complete.\n")
  
  # --- 2. Prepare mRNA Matrix ---
  mrna_mat <- as.data.frame(mrna_df)
  rownames(mrna_mat) <- mrna_mat[[1]]
  mrna_mat[[1]] <- NULL
  colnames(mrna_mat) <- substr(colnames(mrna_mat), 1, 12)
  mrna_mat <- as.matrix(mrna_mat)
  
  # --- 3. Filter and Find Top 500 Genes ---
  mrna_patients_to_keep <- intersect(colnames(mrna_mat), all_patient_ids)
  cat(paste("Found", length(mrna_patients_to_keep), "of your", length(all_patient_ids), "patients in the mRNA data.\n"))
  
  mrna_mat_filtered <- mrna_mat[, mrna_patients_to_keep]
  mrna_vars <- apply(mrna_mat_filtered, 1, var)
  top_mrna <- names(sort(mrna_vars, decreasing = TRUE))[1:500]
  
  # ** This is the unscaled matrix pheatmap will use **
  mrna_mat_top <- mrna_mat_filtered[top_mrna, ]
  
  # --- 4. Scale Data for PAM Clustering ---
  # PAM (like k-means) needs scaled data where patients are rows.
  mrna_scaled <- t(scale(t(mrna_mat_top)))
  mrna_scaled[is.na(mrna_scaled)] <- 0
  
  # ** This is the scaled matrix for clustering **
  mrna_data_for_pam <- t(mrna_scaled) 
  
  # --- 5. Run Iterative PAM Clustering and Plotting ---
  
  # Define a diverging color palette (good for Z-scores from scale="row")
  my_div_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  set.seed(123) # For reproducible PAM clustering
  
  # Loop from k=2 to k=6
  for (k_value in 2:6) {
    
    cat(paste("--- Generating heatmap for k =", k_value, "---\n"))
    
    # 1. Cluster the patients into k groups
    pam_fit <- pam(mrna_data_for_pam, k = k_value)
    
    # 2. Create the column annotation dataframe
    annotation_col <- data.frame(
      Cluster = as.factor(pam_fit$clustering)
    )
    rownames(annotation_col) <- colnames(mrna_mat_top)
    
    # 3. Sort the annotation by cluster number
    annotation_col_sorted <- annotation_col %>%
      arrange(Cluster)
    
    # 4. Get the patient IDs in the new sorted order
    patient_order <- rownames(annotation_col_sorted)
    
    # 5. Generate the heatmap
    # 
    pheatmap(
      mrna_mat_top[, patient_order],  # <-- Use UN-scaled data, sorted by cluster
      color = my_div_colors,        # Use diverging red/blue colors
      scale = "row",                # <-- IMPORTANT: pheatmap scales the rows (Z-score)
      cluster_rows = TRUE,          # Cluster the genes (rows)
      cluster_cols = FALSE,         # DO NOT cluster columns (we already sorted them)
      show_colnames = FALSE,        # Hide patient IDs
      show_rownames = FALSE,        # Hide gene IDs
      main = paste("mRNA Heatmap: Top 500 Genes (PAM k =", k_value, "Clusters)"),
      annotation_col = annotation_col_sorted,  # Add the cluster annotation bar
      gaps_col = find_gaps(annotation_col_sorted, "Cluster") # Add vertical lines
    )
  }
  
  
  # --- 1. Load Libraries ---
  # (We assume pheatmap and clValid are already installed from before)
  library(pheatmap)
  library(clValid)
  
  
  # --- 2. Download Copy Number Data (Direct Method) ---
  cat("Downloading TCGA-LIHC.gene-level_ascat2.tsv.gz...\n")
  url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.gene-level_ascat2.tsv.gz"
  dest_file <- "TCGA-LIHC.gene-level_ascat2.tsv.gz"
  
  # Download the file from the URL
  download.file(url, dest_file)
  
  # Read the compressed file directly into R
  cnv_data <- read.delim(gzfile(dest_file), 
                         row.names = 1, 
                         check.names = FALSE)
  
  cat("Download complete.\n")
  # print(dim(cnv_data)) # Should be 60624 x 372
  
  
  # --- 3. Data Preparation ---
  
  ### Find Top 500 Variable Genes
  # Calculate variance for each gene (row)
  row_vars <- apply(cnv_data, 1, var)
  
  # Filter out any genes with 0 variance (which can cause errors)
  row_vars <- row_vars[!is.na(row_vars) & row_vars > 0]
  
  # Get the names of the top 500 most variable genes
  top_500_genes <- names(sort(row_vars, decreasing = TRUE)[1:500])
  
  # Create a new matrix with only these 500 genes
  data_subset <- cnv_data[top_500_genes, ]
  
  ### Row-Normalize Data (Z-score)
  # We scale the data for each gene (row) to have a mean of 0 and SD of 1
  # This highlights relative amplifications/deletions for each gene
  scaled_data <- t(scale(t(data_subset)))
  
  # Handle any potential NaN values
  scaled_data[is.na(scaled_data)] <- 0
  
  
  # --- 4. Heatmap Visualization (k=2 and k=3) ---
  
  ### Heatmap for k=2 Clusters
  cat("Generating heatmap for k=2... (this may take a moment)\n")
  pheatmap(
    scaled_data,
    main = "Top 500 Variable Genes (Copy Number, k=2)",
    scale = "none",          # Data is already scaled
    show_colnames = FALSE,   # Hides the sample names
    show_rownames = FALSE,   # Hides the gene names
    cluster_cols = TRUE,     # Cluster the samples
    cluster_rows = TRUE,     # Cluster the genes
    cutree_cols = 2          # Cut sample cluster tree into 2 groups
  )
  
  
 
  
  ### Heatmap for k=3 Clusters
  cat("Generating heatmap for k=3...\n")
  pheatmap(
    scaled_data,
    main = "Top 500 Variable Genes (Copy Number, k=3)",
    scale = "none",
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    cutree_cols = 3          # Cut sample cluster tree into 3 groups
  )
  
  
  # --- 5. Internal Cluster Validation (k=2 to 6) ---
  # We use k-means on the transposed (samples in rows) data
  
  cat("Running internal cluster validation for k=2 to 6...\n")
  # This step may take a minute
  cl_validation <- clValid(
    obj = t(scaled_data),         # Data transposed (samples in rows)
    nClust = 2:6,                 # k values to test
    clMethods = "kmeans",         # K-Means algorithm
    validation = "internal"
  )
  
  # Print the results
  cat("--- Cluster Validation Results (Copy Number) ---\n")
  print(summary(cl_validation))
  
  cat("\n--- Optimal Scores (Copy Number) ---\n")
  print(optimalScores(cl_validation))
  
  
  # --- 6. Plot the Validation Metrics ---
  cat("Generating cluster validation plots...\n")
  
  # Save the plot to a new file
  png("cnv_cluster_validation_plots.png", width = 800, height = 800)
  plot(cl_validation)
  dev.off()
  
  cat("Saved 'cnv_cluster_validation_plots.png' to your working directory.\n")
  
  
  
  
  
  
  
  
  
  
  library(data.table)
  
  
  # --- 2. Read the file with fread (The Fast Way) ---
  cat("Reading data file with data.table::fread()...\n")
  dest_file <- "TCGA-LIHC.methylation450.tsv.gz"
  
  # 'data.table = FALSE' tells it to return a regular data.frame
  # 'check.names = FALSE' keeps sample IDs like "TCGA-..." intact
  meth_data_df <- fread(dest_file, 
                        data.table = FALSE, 
                        check.names = FALSE)
  
  # --- 3. Set Row Names ---
  # fread() reads the probe IDs as the first column.
  # We need to manually set them as row names and then remove the first column.
  rownames(meth_data_df) <- meth_data_df[, 1]
  meth_data <- meth_data_df[, -1]
  
  cat("Read complete.\n")
  
  print(head(meth_data))
  
  
  cat("Calculating variance for 486,000+ probes...\n")
  row_vars <- apply(meth_data, 1, var) 
  
  # Filter out any probes with 0 variance (which can cause errors)
  row_vars <- row_vars[!is.na(row_vars) & row_vars > 0]
  
  # Get the names of the top 500 most variable probes
  top_500_probes <- names(sort(row_vars, decreasing = TRUE)[1:500])
  
  # Create a new matrix with only these 500 probes
  data_subset <- meth_data[top_500_probes, ]
  cat("Top 500 probes selected.\n")
  
  ### Row-Normalize Data (Z-score)
  # We scale the data for each probe (row) to have a mean of 0 and SD of 1
  scaled_data <- t(scale(t(data_subset)))
  
  # Handle any potential NaN values
  scaled_data[is.na(scaled_data)] <- 0
  
  
  # --- 4. Heatmap Visualization (k=2 and k=3) ---
  
  ### Heatmap for k=2 Clusters
  cat("Generating heatmap for k=2...\n")
  pheatmap(
    scaled_data,
    main = "Top 500 Variable Probes (Methylation, k=2)",
    scale = "none",          # Data is already scaled
    show_colnames = FALSE,   # Hides the sample names
    show_rownames = FALSE,   # Hides the probe names
    cluster_cols = TRUE,     # Cluster the samples
    cluster_rows = TRUE,     # Cluster the genes
    cutree_cols = 2          # Cut sample cluster tree into 2 groups
  )
  
  

  
  
  ### Heatmap for k=3 Clusters
  cat("Generating heatmap for k=3...\n")
  pheatmap(
    scaled_data,
    main = "Top 500 Variable Probes (Methylation, k=3)",
    scale = "none",
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    cutree_cols = 3          # Cut sample cluster tree into 3 groups
  )
  
  
  # --- 5. Internal Cluster Validation (k=2 to 6) ---
  # We use k-means on the transposed (samples in rows) data
  library(clValid)
  cat("Running internal cluster validation for k=2 to 6...\n")
  cl_validation <- clValid(
    obj = t(scaled_data),         # Data transposed (samples in rows)
    nClust = 2:6,                 # k values to test
    clMethods = "kmeans",         # K-Means algorithm
    validation = "internal"
  )
  
  # Print the results
  cat("--- Cluster Validation Results (Methylation) ---\n")
  print(summary(cl_validation))
  
  cat("\n--- Optimal Scores (Methylation) ---\n")
  print(optimalScores(cl_validation))
  
  
  # --- 6. Plot the Validation Metrics ---
  cat("Generating cluster validation plots...\n")
  
  # Save the plot to a new file
  png("methylation_cluster_validation_plots.png", width = 800, height = 800)
  plot(cl_validation)
  dev.off()
  
  cat("Saved 'methylation_cluster_validation_plots.png' to your working directory.\n")
  
  