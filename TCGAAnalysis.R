# First, ensure you have the BiocManager package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# --- 2. Install TCGAbiolinks, if you don't have it ---
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require(UCSCXenaTools)) install.packages("UCSCXenaTools")
library(UCSCXenaTools)
# Now, use BiocManager to install XenaR
BiocManager::install("XenaR")

install.packages(c(
  "dplyr",        # Provides the %>% pipe operator and data wrangling
  "data.table",   # Provides the fast fread() function
  "stringr",      # For stringr::str_extract() (Cleaning Method 3)
  "pheatmap",     # For generating the heatmaps
  "factoextra",   # For fviz_nbclust() (Elbow/Silhouette plots)
  "ggplot2",      # For plotting (ggsave)
  "RColorBrewer"  # For the heatmap color palettes
))
library(dplyr)
library(data.table)

# --- 3. Load the library ---
# This is the line that fixes your error
library(TCGAbiolinks)
library(DESeq2)
library(pheatmap)
library(ConsensusClusterPlus)
library(RColorBrewer)
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




#--------------------------------mrna------------------------



dataset_mrna <- "TCGA.LIHC.sampleMap/HiSeqV2_PANCAN"
hub_mrna <- "https://tcga.xenahubs.net" # Using the hub you specified

cat("--- Downloading mRNA data... ---\n")

# ------------------------------------------------------------
# 3. Query and Download the mRNA data
# ------------------------------------------------------------
xe_mrna <- XenaGenerate() %>%
  XenaFilter(filterDatasets = dataset_mrna) %>%
  XenaFilter(filterHosts = hub_mrna) %>%
  XenaQuery() %>%
  XenaDownload(download_R = FALSE) # Set to FALSE, we will read with fread

# ------------------------------------------------------------
# 4. Load mRNA matrix into R
# ------------------------------------------------------------
# Find the downloaded file path
filepath_mrna <- list.files(
  path = tempdir(),
  pattern = "HiSeqV2_PANCAN\\.gz$", # Pattern from the dataset ID
  full.names = TRUE,
  recursive = TRUE
)

cat(paste("Reading file:", filepath_mrna[1], "\n"))

# Read the matrix using data.table for speed
mrna_df <- fread(filepath_mrna[1], sep = "\t", header = TRUE)

# ------------------------------------------------------------
# 5. Wrangle the mRNA matrix
# ------------------------------------------------------------
# This dataset has 423 samples, including normals (e.g., -11)
# We must clean the column names to 12-char patient IDs to match your lists

# Get all column names except the first (which is 'sample')
mrna_cols <- colnames(mrna_df)
mrna_cols_clean <- mrna_cols

# Apply the 12-character substring to all sample columns
mrna_cols_clean[-1] <- substr(mrna_cols[-1], 1, 12)

# Set the cleaned names back
colnames(mrna_df) <- mrna_cols_clean
print(rownames(mrna_mat))
# Convert to a numeric matrix with gene symbols as rownames
# (This follows the same logic as your miRNA code)
mrna_mat <- as.data.frame(mrna_df)
rownames(mrna_mat) <- mrna_mat[[1]] 
print(rownames(mrna_mat))# Set genes (e.g., "?|100130426") as rownames
mrna_mat[[1]] <- NULL              # Remove the gene ID column
mrna_mat <- as.matrix(mrna_mat)    # Convert to a numeric matrix

cat("--- Original mRNA Matrix Dimensions ---\n")
print(dim(mrna_mat)) # Should be 20531 genes x 423 samples

# ------------------------------------------------------------
# 6. Filter matrix to *only* your 373 patients
# ------------------------------------------------------------
# 'all_patient_ids' should exist from your previous code block

# Find the intersection of patients in the mRNA data and your cluster list
patients_to_keep_mrna <- intersect(colnames(mrna_mat), all_patient_ids)

# Print a summary of the overlap
cat(paste("\nFound", length(patients_to_keep_mrna), "of your", 
          length(all_patient_ids), "patients in the mRNA data.\n"))

# Subset the matrix to *only* these patients (columns)
mrna_mat_filtered <- mrna_mat[, patients_to_keep_mrna]

# Inspect new dimensions
cat("--- Filtered mRNA Matrix Dimensions ---\n")
print(dim(mrna_mat_filtered)) # Should be 20531 genes x 373 patients

cat("\n--- Head of filtered mRNA matrix ---\n")
print(mrna_mat_filtered[1:4, 1:4])



# ------------------------------------------------------------
# 1. Install/Load Libraries
# ------------------------------------------------------------
# We need pheatmap for the heatmaps
if (!require(pheatmap)) install.packages("pheatmap")
# We need factoextra to calculate optimal k scores
if (!require(factoextra)) install.packages("factoextra")
# For the color palettes
if (!require(RColorBrewer)) install.packages("RColorBrewer")

library(pheatmap)
library(factoextra)
library(RColorBrewer)

# ------------------------------------------------------------
# 2. Prepare mRNA data for clustering
# ------------------------------------------------------------
# 'mrna_mat_filtered' should exist from the previous code block
# It has 20531 genes (rows) x 373 patients (columns)

cat("--- Preparing mRNA data for clustering ---\n")

# --- Step 2a: Select top 500 most variable genes ---
# Calculate variance for each gene (row-wise)
gene_vars <- apply(mrna_mat_filtered, 1, var)

# Get the names of the top 500 genes
top_500_genes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# Subset the matrix
mrna_mat_top <- mrna_mat_filtered[top_500_genes, ]
cat(paste("Selected top 500 variable genes. Matrix dimensions:", 
          dim(mrna_mat_top)[1], "genes x", dim(mrna_mat_top)[2], "patients\n"))

# --- Step 2b: Z-score the data (Scale by gene/row) ---
# For k-means and heatmaps, we scale data so each gene (row) has a mean of 0 
# and a standard deviation of 1.
# We use t(scale(t(...))) to scale by row.
mrna_mat_top_scaled <- t(scale(t(mrna_mat_top)))

# Verify scaling
cat("Verifying Z-score (row means should be ~0):\n")
print(head(apply(mrna_mat_top_scaled, 1, mean)))

# --- Step 2c: Transpose for k-means ---
# k-means clusters rows, so we need patients as rows and genes as columns
data_for_kmeans_mrna <- t(mrna_mat_top_scaled)
cat(paste("\nTransposed matrix for k-means. Dimensions:",
          dim(data_for_kmeans_mrna)[1], "patients x", 
          dim(data_for_kmeans_mrna)[2], "genes\n"))


# ------------------------------------------------------------
# 3. Calculate Optimal 'k' (Scores)
# ------------------------------------------------------------
# Set seed for reproducible results
set.seed(123)

cat("--- Calculating Optimal 'k' (k=2 to 6) ---\n")

# --- Elbow Method (WSS) ---
cat("Generating Elbow method plot...\n")

# 1. Create the plot object
fviz_nbclust(
  data_for_kmeans_mrna, 
  kmeans, 
  method = "wss", # "wss" = Within Sum of Squares
  k.max = 6,
  nstart = 25
)


# --- Silhouette Method ---
cat("Generating Silhouette method plot...\n")
# This creates the plot object AND prints it in one step
print(head(data_for_kmeans_mrna))
# 1. Create the plot object
 fviz_nbclust(
  data_for_kmeans_mrna,
  kmeans,
  method = "silhouette",
  k.max = 6,
  nstart = 25
) 

# Optional: still print to RStudio screen

cat("Optimal 'k' plots saved.\n")

# ------------------------------------------------------------
# 4. Run K-Means for k=2 and k=3
# ------------------------------------------------------------
set.seed(123)
cat("--- Running K-Means for k=2 and k=3 ---\n")
kmeans_k2_mrna <- kmeans(data_for_kmeans_mrna, centers = 2, nstart = 25)
kmeans_k3_mrna <- kmeans(data_for_kmeans_mrna, centers = 3, nstart = 25)
cat("K-Means clustering complete.\n")

# ------------------------------------------------------------
# 5. Print Patient IDs for k=3
# ------------------------------------------------------------
cat("\n\n--- PATIENT IDs FOR K=3 CLUSTERS ---\n")
cluster_assignments_k3 <- kmeans_k3_mrna$cluster
all_patient_names <- names(cluster_assignments_k3)

for (i in 1:3) {
  cat(paste("\n--- mRNA K-Means Cluster", i, "(k=3) ---"))
  patients_in_cluster <- all_patient_names[cluster_assignments_k3 == i]
  cat(paste("\nTotal patients:", length(patients_in_cluster), "\n"))
  cat(patients_in_cluster, sep = "\n")
}

cluster1_ids <- all_patient_names[cluster_assignments_k3 == 1]
cluster2_ids <- all_patient_names[cluster_assignments_k3 == 2]
cluster3_ids <- all_patient_names[cluster_assignments_k3 == 3]

all_patient_ids <- c(cluster1_ids, 
                     cluster2_ids, 
                     cluster3_ids)

# ... (Verification print statements are fine) ...

# ------------------------------------------------------------
# 6. Create Annotation Data Frame
# ------------------------------------------------------------
# This data frame looks correct with your new 'clustermrna' column
annotation_df_mrna <- data.frame(
  KMeans_k2 = as.factor(kmeans_k2_mrna$cluster),
  KMeans_k3 = as.factor(kmeans_k3_mrna$cluster),
  clustermrna = as.factor(kmeans_k3_mrna$cluster),
  row.names = all_patient_names
)
cat("\n--- Annotation Data Frame for Heatmaps (Head) ---\n")
print(head(annotation_df_mrna))

# ------------------------------------------------------------
# 7. Generate Heatmaps
# ------------------------------------------------------------
my_heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# --- Heatmap 1: K=2 ---
cat("\nGenerating k=2 heatmap (saving to 'kmeans_mrna_k2_heatmap.png')...\n")
col_order_k2 <- order(annotation_df_mrna$KMeans_k2)
mat_ordered_k2 <- mrna_mat_top_scaled[, col_order_k2]

# --- ADJUSTMENT ---
# Select ALL columns from the annotation frame, not just one
anno_ordered_k2 <- annotation_df_mrna[col_order_k2, , drop = FALSE]

pheatmap(
  mat_ordered_k2,
  scale = "none",
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k2, # Pass the multi-column annotation
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC mRNA: K-Means (k=2) Clusters (Top 500 Var Genes)",
  width = 12,
  height = 10
)


# --- Heatmap 2: K=3 ---
cat("Generating k=3 heatmap (saving to 'kmeans_mrna_k3_heatmap.png')...\n")

# --- ADJUSTMENT ---
# We use 'all_kmeans_patient_ids' from your previous code chunk.
# This list is ALREADY sorted by Cluster 1, then 2, then 3.
# This assumes 'all_kmeans_patient_ids' exists in your environment.
if (!exists("all_kmeans_patient_ids")) {
  stop("Error: 'all_kmeans_patient_ids' not found. Please run the previous code chunk.")
}
col_order_k3 <- all_patient_ids

# Order the matrix and annotation using this pre-sorted list
mat_ordered_k3 <- mrna_mat_top_scaled[, col_order_k3]

# --- ADJUSTMENT ---
# Select ALL annotation columns, ordered by our new 'col_order_k3'
anno_ordered_k3 <- annotation_df_mrna[col_order_k3, , drop = FALSE]

pheatmap(
  mat_ordered_k3,
  scale = "none",
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k3, # Pass the multi-column annotation
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC mRNA: K-Means (k=3) Clusters (Top 500 Var Genes)",
  width = 12,
  height = 10
)

cat("\n--- All tasks complete. Check your working directory for PNG files. ---\n")








dataset_mirna <- "TCGA.LIHC.sampleMap/miRNA_HiSeq_gene"
hub_mirna <- "https://tcga.xenahubs.net" 

cat("--- Downloading miRNA data... ---\n")

xe_mirna <- XenaGenerate() %>%
  XenaFilter(filterDatasets = dataset_mirna) %>%
  XenaFilter(filterHosts = hub_mirna) %>%
  XenaQuery() %>%
  XenaDownload(download_R = FALSE) # Set to FALSE, we will read with fread

# ------------------------------------------------------------
# 4. Load miRNA matrix into R
# ------------------------------------------------------------
# Find the downloaded file path
filepath_mirna <- list.files(
  path = tempdir(),
  pattern = "miRNA_HiSeq_gene\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)

cat(paste("Reading file:", filepath_mirna[1], "\n"))
mirna_df <- fread(filepath_mirna[1], sep = "\t", header = TRUE)

# ------------------------------------------------------------
# 5. Wrangle the miRNA matrix
# ------------------------------------------------------------
mirna_cols <- colnames(mirna_df)
mirna_cols_clean <- mirna_cols
mirna_cols_clean[-1] <- substr(mirna_cols[-1], 1, 12) # Clean column names
colnames(mirna_df) <- mirna_cols_clean

# Convert to a numeric matrix
mirna_mat <- as.data.frame(mirna_df)
rownames(mirna_mat) <- mirna_mat[[1]] 
mirna_mat[[1]] <- NULL               
mirna_mat <- as.matrix(mirna_mat)     

cat("--- Original miRNA Matrix Dimensions ---\n")
print(dim(mirna_mat)) 

# ------------------------------------------------------------
# 6. Filter matrix to *only* your 373 patients
# ------------------------------------------------------------
# 'all_patient_ids' should exist from your previous script
patients_to_keep_mirna <- intersect(colnames(mirna_mat), all_patient_ids)

cat(paste("\nFound", length(patients_to_keep_mirna), "of your", 
          length(all_patient_ids), "patients in the miRNA data.\n"))

mirna_mat_filtered <- mirna_mat[, patients_to_keep_mirna]

cat("--- Filtered miRNA Matrix Dimensions ---\n")
print(dim(mirna_mat_filtered))
cat("\n--- Head of filtered miRNA matrix ---\n")
print(mirna_mat_filtered[1:4, 1:4])


# ------------------------------------------------------------
# 7. Prepare miRNA data for clustering (Robust Version)
# ------------------------------------------------------------
cat("--- Preparing miRNA data for clustering ---\n")

# --- Step 7a: Select top 500 most variable miRNAs ---
# (This step was from the previous fix, it's still correct)
cat("Calculating variance...\n")
mirna_vars <- apply(mirna_mat_filtered, 1, var, na.rm = TRUE)
mirna_vars_clean <- mirna_vars[!is.na(mirna_vars)]
mirna_vars_sorted <- sort(mirna_vars_clean, decreasing = TRUE)
top_500_mirnas <- names(head(mirna_vars_sorted, 500))
mirna_mat_top <- mirna_mat_filtered[top_500_mirnas, ]

cat(paste("Selected top", length(top_500_mirnas), "variable miRNAs.\n"))

# --- Step 7b: Z-score the data (Scale by miRNA/row) ---
# We must scale manually to correctly handle NA values
cat("Scaling data (Z-score by row) and ignoring NAs...\n")
row_means <- rowMeans(mirna_mat_top, na.rm = TRUE)
row_sds <- apply(mirna_mat_top, 1, sd, na.rm = TRUE)
mirna_mat_top_scaled <- (mirna_mat_top - row_means) / row_sds

# --- Step 7c: Remove problematic miRNAs (NA/NaN/Inf) ---
# Rows with sd=0 will now be Inf or NaN. We also check for any remaining NAs.
# !is.finite() checks for NA, NaN, Inf, and -Inf all at once.
rows_with_nonfinite <- apply(mirna_mat_top_scaled, 1, function(row) any(!is.finite(row)))

if(any(rows_with_nonfinite)) {
  cat(paste("Found and removed", sum(rows_with_nonfinite), "miRNA(s) with NA, NaN, or Inf values.\n"))
  # Keep only the "clean" rows (miRNAs)
  mirna_mat_top_scaled_clean <- mirna_mat_top_scaled[!rows_with_nonfinite, ]
} else {
  cat("No problematic miRNAs (NA/NaN/Inf) found.\n")
  mirna_mat_top_scaled_clean <- mirna_mat_top_scaled
}

# --- Step 7d: Transpose for k-means ---
data_for_kmeans_mirna_raw <- t(mirna_mat_top_scaled_clean)

# --- Step 7e (FINAL CHECK): Remove problematic patients ---
# k-means cannot handle *any* NA/NaN/Inf. It's possible a patient
# had NA for all the miRNAs we kept. We must remove them.
patients_with_nonfinite <- apply(data_for_kmeans_mirna_raw, 1, function(row) any(!is.finite(row)))

if(any(patients_with_nonfinite)) {
  cat(paste("Found and removed", sum(patients_with_nonfinite), "patient(s) with NA/NaN/Inf values.\n"))
  # This is the final, clean data for kmeans
  data_for_kmeans_mirna <- data_for_kmeans_mirna_raw[!patients_with_nonfinite, ]
} else {
  cat("All patients have finite data. Ready for clustering.\n")
  data_for_kmeans_mirna <- data_for_kmeans_mirna_raw
}

cat(paste("\nTransposed matrix for k-means. Final dimensions:",
          dim(data_for_kmeans_mirna)[1], "patients x", 
          dim(data_for_kmeans_mirna)[2], "miRNAs\n"))


# ------------------------------------------------------------
# 8. Calculate Optimal 'k' (Scores)
# ------------------------------------------------------------
set.seed(123)
cat("--- Calculating Optimal 'k' (k=2 to 6) for miRNA ---\n")

# --- Elbow Method (WSS) ---
cat("Generating Elbow method plot...\n")
# (FIXED: Use the correct data object)
p_elbow_mirna <- fviz_nbclust(
  data_for_kmeans_mirna, 
  kmeans, 
  method = "wss",
  k.max = 6,
  nstart = 25
) + ggtitle("miRNA K-Means: Elbow Method")
print(p_elbow_mirna)


# --- Silhouette Method ---
cat("Generating Silhouette method plot...\n")
# (FIXED: Use the correct data object)
p_sil_mirna <- fviz_nbclust(
  data_for_kmeans_mirna,
  kmeans,
  method = "silhouette",
  k.max = 6,
  nstart = 25
) + ggtitle("miRNA K-Means: Silhouette Method")
print(p_sil_mirna)

cat("Optimal 'k' plots displayed.\n")

# ------------------------------------------------------------
# 9. Run K-Means for k=2 and k=3
# ------------------------------------------------------------
set.seed(123)
cat("--- Running K-Means for k=2 and k=3 ---\n")
# (FIXED: Use the correct data object)
kmeans_k2_mirna <- kmeans(data_for_kmeans_mirna, centers = 2, nstart = 25)
kmeans_k3_mirna <- kmeans(data_for_kmeans_mirna, centers = 3, nstart = 25)
cat("K-Means clustering complete.\n")

# ------------------------------------------------------------
# 10. Print Patient IDs for k=3
# ------------------------------------------------------------
cat("\n\n--- PATIENT IDs FOR K=3 CLUSTERS (miRNA) ---\n")
cluster_assignments_k3_mirna <- kmeans_k3_mirna$cluster
all_patient_names_mirna <- names(cluster_assignments_k3_mirna)

for (i in 1:3) {
  cat(paste("\n--- miRNA K-Means Cluster", i, "(k=3) ---"))
  patients_in_cluster <- all_patient_names_mirna[cluster_assignments_k3_mirna == i]
  cat(paste("\nTotal patients:", length(patients_in_cluster), "\n"))
  cat(patients_in_cluster, sep = "\n")
}

kmeans_cluster1_ids_mirna <- all_patient_names_mirna[cluster_assignments_k3_mirna == 1]
kmeans_cluster2_ids_mirna <- all_patient_names_mirna[cluster_assignments_k3_mirna == 2]
kmeans_cluster3_ids_mirna <- all_patient_names_mirna[cluster_assignments_k3_mirna == 3]

all_kmeans_patient_ids_mirna <- c(kmeans_cluster1_ids_mirna, 
                                  kmeans_cluster2_ids_mirna, 
                                  kmeans_cluster3_ids_mirna)

# ------------------------------------------------------------
# 11. Create Annotation Data Frame
# ------------------------------------------------------------
# ------------------------------------------------------------
# 6. Create Annotation Data Frames
# ------------------------------------------------------------
# You may need to load dplyr if you haven't already
# library(dplyr) 

# --- Step 6a: Create the miRNA annotation data ---
# (This has the patients from the miRNA k-means)
mirna_anno_df <- data.frame(
  patient_id = all_patient_names_mirna, # Use patient_id as a COLUMN
  KMeans_miRNA_k2 = as.factor(kmeans_k2_mirna$cluster),
  KMeans_miRNA_k3 = as.factor(kmeans_k3_mirna$cluster),
  clustermiRNA = as.factor(kmeans_k3_mirna$cluster) # You can keep this
)

# --- Step 6b: Create the mRNA annotation data ---
# (This has the patients from the mRNA k-means)
mrna_anno_df <- data.frame(
  patient_id = names(kmeans_k3_mrna$cluster), # Get IDs from the cluster object
  clustermrna = as.factor(kmeans_k3_mrna$cluster)
  # You could also add the k=2 mRNA clusters here if you want
  # KMeans_mRNA_k2 = as.factor(kmeans_k2_mrna$cluster)
)

# --- Step 6c: Join them using the patient_id "glue" ---
# A "left_join" keeps everything from the "left" data (miRNA)
# and adds matching data from the "right" data (mRNA).
# Any miRNA patient not in the mRNA analysis will get 'NA'
annotation_df_mirna <- left_join(
  mirna_anno_df, 
  mrna_anno_df, 
  by = "patient_id" # The common column to match on
)

# --- Step 6d: Set row names (required for pheatmap) ---
# Now that it's merged, set the row names from the patient_id column
rownames(annotation_df_mirna) <- annotation_df_mirna$patient_id

cat("\n--- Merged Annotation Data Frame for Heatmaps (Head) ---\n")
print(head(annotation_df_mirna))

# ------------------------------------------------------------
# 12. Generate Heatmaps
# ------------------------------------------------------------
my_heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# --- Heatmap 1: K=2 ---
cat("\nGenerating k=2 heatmap (miRNA)...\n")
col_order_k2 <- order(annotation_df_mirna$KMeans_miRNA_k2)
# (FIXED: Use the clean, scaled matrix)
mat_ordered_k2 <- mirna_mat_top_scaled_clean[, col_order_k2]
anno_ordered_k2 <- annotation_df_mirna[col_order_k2, , drop = FALSE]

pheatmap(
  mat_ordered_k2,
  scale = "none",
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC miRNA: K-Means (k=2) Clusters (Top 500 Var miRNAs)"
)

# --- Heatmap 2: K=3 (Ordered by Cluster 1, 2, 3) ---
cat("Generating k=3 heatmap (miRNA)...\n")
col_order_k3 <- all_kmeans_patient_ids_mirna
# (FIXED: Use the clean, scaled matrix)
mat_ordered_k3 <- mirna_mat_top_scaled_clean[, col_order_k3]
anno_ordered_k3 <- annotation_df_mirna[col_order_k3, , drop = FALSE]

pheatmap(
  mat_ordered_k3,
  scale = "none",
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k3,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC miRNA: K-Means (k=3) Clusters (Top 500 Var miRNAs)"
)

cat("\n--- All tasks complete. Check your R Plots pane. ---\n")










dataset_cnv <- "TCGA.LIHC.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
hub_cnv <- "https://tcga.xenahubs.net" 

cat("--- Downloading CNV data... ---\n")

xe_cnv <- XenaGenerate() %>%
  XenaFilter(filterDatasets = dataset_cnv) %>%
  XenaFilter(filterHosts = hub_cnv) %>%
  XenaQuery() %>%
  XenaDownload(download_R = FALSE) # Set to FALSE, we will read with fread

# ------------------------------------------------------------
# 4. Load CNV matrix into R
# ------------------------------------------------------------
# Find the downloaded file path
filepath_cnv <- list.files(
  path = tempdir(),
  pattern = "Gistic2_CopyNumber_Gistic2_all_thresholded\\.by_genes\\.gz$", # <-- CHANGED PATTERN
  full.names = TRUE,
  recursive = TRUE
)

cat(paste("Reading file:", filepath_cnv[1], "\n"))

# Read the matrix using data.table for speed
cnv_df <- fread(filepath_cnv[1], sep = "\t", header = TRUE)

# ------------------------------------------------------------
# 5. Wrangle the CNV matrix
# ------------------------------------------------------------
# This dataset has 370 samples
cnv_cols <- colnames(cnv_df)
cnv_cols_clean <- cnv_cols
cnv_cols_clean[-1] <- substr(cnv_cols[-1], 1, 12) # Clean column names
colnames(cnv_df) <- cnv_cols_clean

# Convert to a numeric matrix with gene symbols as rownames
cnv_mat <- as.data.frame(cnv_df)
rownames(cnv_mat) <- cnv_mat[[1]] # Set gene IDs as rownames
cnv_mat[[1]] <- NULL             # Remove the ID column
cnv_mat <- as.matrix(cnv_mat)     # Convert to a numeric matrix

cat("--- Original CNV Matrix Dimensions ---\n")
print(dim(cnv_mat)) # Should be 24777 genes x 370 samples

# ------------------------------------------------------------
# 6. Filter matrix to *only* your 373 patients
# ------------------------------------------------------------
patients_to_keep_cnv <- intersect(colnames(cnv_mat), all_patient_ids)

cat(paste("\nFound", length(patients_to_keep_cnv), "of your", 
          length(all_patient_ids), "patients in the CNV data.\n"))

# Subset the matrix
cnv_mat_filtered <- cnv_mat[, patients_to_keep_cnv]

cat("--- Filtered CNV Matrix Dimensions ---\n")
print(dim(cnv_mat_filtered)) # ~24777 genes x 370 patients

cat("\n--- Head of filtered CNV matrix ---\n")
print(cnv_mat_filtered[1:4, 1:4])


# ------------------------------------------------------------
# 7. Prepare CNV data for clustering (Robust Version)
# ------------------------------------------------------------
# We follow the same robust process as before
cat("--- Preparing CNV data for clustering ---\n")

# --- Step 7a: Select top 500 most variable genes ---
cat("Calculating variance...\n")
cnv_vars <- apply(cnv_mat_filtered, 1, var, na.rm = TRUE)
cnv_vars_clean <- cnv_vars[!is.na(cnv_vars)]
cnv_vars_sorted <- sort(cnv_vars_clean, decreasing = TRUE)
top_500_genes_cnv <- names(head(cnv_vars_sorted, 500))
cnv_mat_top <- cnv_mat_filtered[top_500_genes_cnv, ]

cat(paste("Selected top", length(top_500_genes_cnv), "variable genes.\n"))

# --- Step 7b: Z-score the data (Scale by gene/row) ---
cat("Scaling data (Z-score by row) and ignoring NAs...\n")
row_means <- rowMeans(cnv_mat_top, na.rm = TRUE)
row_sds <- apply(cnv_mat_top, 1, sd, na.rm = TRUE)
cnv_mat_top_scaled <- (cnv_mat_top - row_means) / row_sds

# --- Step 7c: Remove problematic genes (NA/NaN/Inf) ---
# !is.finite() checks for NA, NaN, Inf
rows_with_nonfinite <- apply(cnv_mat_top_scaled, 1, function(row) any(!is.finite(row)))

if(any(rows_with_nonfinite)) {
  cat(paste("Found and removed", sum(rows_with_nonfinite), "gene(s) with NA, NaN, or Inf values (likely 0 variance).\n"))
  cnv_mat_top_scaled_clean <- cnv_mat_top_scaled[!rows_with_nonfinite, ]
} else {
  cat("No problematic genes (NA/NaN/Inf) found.\n")
  cnv_mat_top_scaled_clean <- cnv_mat_top_scaled
}

# --- Step 7d: Transpose for k-means ---
data_for_kmeans_cnv_raw <- t(cnv_mat_top_scaled_clean)

# --- Step 7e (FINAL CHECK): Remove problematic patients ---
patients_with_nonfinite <- apply(data_for_kmeans_cnv_raw, 1, function(row) any(!is.finite(row)))

if(any(patients_with_nonfinite)) {
  cat(paste("Found and removed", sum(patients_with_nonfinite), "patient(s) with NA/NaN/Inf values.\n"))
  data_for_kmeans_cnv <- data_for_kmeans_cnv_raw[!patients_with_nonfinite, ]
} else {
  cat("All patients have finite data. Ready for clustering.\n")
  data_for_kmeans_cnv <- data_for_kmeans_cnv_raw
}

cat(paste("\nTransposed matrix for k-means. Final dimensions:",
          dim(data_for_kmeans_cnv)[1], "patients x", 
          dim(data_for_kmeans_cnv)[2], "genes\n"))


# ------------------------------------------------------------
# 8. Calculate Optimal 'k' (Scores)
# ------------------------------------------------------------
set.seed(123)
cat("--- Calculating Optimal 'k' (k=2 to 6) for CNV ---\n")

# --- Elbow Method (WSS) ---
cat("Generating Elbow method plot...\n")
p_elbow_cnv <- fviz_nbclust(
  data_for_kmeans_cnv, 
  kmeans, 
  method = "wss",
  k.max = 6,
  nstart = 25
) + ggtitle("CNV K-Means: Elbow Method")
print(p_elbow_cnv) # Print to screen


# --- Silhouette Method ---
cat("Generating Silhouette method plot...\n")
p_sil_cnv <- fviz_nbclust(
  data_for_kmeans_cnv,
  kmeans,
  method = "silhouette",
  k.max = 6,
  nstart = 25
) + ggtitle("CNV K-Means: Silhouette Method")
print(p_sil_cnv) # Print to screen

cat("Optimal 'k' plots displayed.\n")

# ------------------------------------------------------------
# 9. Run K-Means for k=2 and k=3
# ------------------------------------------------------------
set.seed(123)
cat("--- Running K-Means for k=2 and k=3 ---\n")
kmeans_k2_cnv <- kmeans(data_for_kmeans_cnv, centers = 2, nstart = 25)
kmeans_k3_cnv <- kmeans(data_for_kmeans_cnv, centers = 3, nstart = 25)
cat("K-Means clustering complete.\n")

# ------------------------------------------------------------
# 10. Get Patient IDs for k=3
# ------------------------------------------------------------
cat("\n\n--- PATIENT IDs FOR K=3 CLUSTERS (CNV) ---\n")
cluster_assignments_k3_cnv <- kmeans_k3_cnv$cluster
all_patient_names_cnv <- names(cluster_assignments_k3_cnv)

# (We just do this to create the sorted list for the k=3 heatmap)
kmeans_cluster1_ids_cnv <- all_patient_names_cnv[cluster_assignments_k3_cnv == 1]
kmeans_cluster2_ids_cnv <- all_patient_names_cnv[cluster_assignments_k3_cnv == 2]
kmeans_cluster3_ids_cnv <- all_patient_names_cnv[cluster_assignments_k3_cnv == 3]

# This is the list sorted by the *newly calculated* k-means clusters
all_kmeans_patient_ids_cnv <- c(kmeans_cluster1_ids_cnv, 
                                kmeans_cluster2_ids_cnv, 
                                kmeans_cluster3_ids_cnv)

# ------------------------------------------------------------
# 11. Create Merged Annotation Data Frame
# ------------------------------------------------------------
# This block assumes 'kmeans_k3_mrna' and 'kmeans_k3_mirna'
# exist in your environment from previous runs.
cat("--- Building Merged Annotation Data Frame ---\n")

# --- Step 11a: Create the NEW CNV annotation data ---
cnv_anno_df <- data.frame(
  patient_id = all_patient_names_cnv, # Use patient_id as a COLUMN
  KMeans_CNV_k2 = as.factor(kmeans_k2_cnv$cluster),
  KMeans_CNV_k3 = as.factor(kmeans_k3_cnv$cluster)
)

# --- Step 11b: Create data frame for OLD mRNA clusters ---
if (exists("kmeans_k3_mrna")) {
  mrna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mrna$cluster),
    clustermrna = as.factor(kmeans_k3_mrna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mrna' object not found. Skipping mRNA annotation.\n")
  mrna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11c: Create data frame for OLD miRNA clusters ---
if (exists("kmeans_k3_mirna")) {
  mirna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mirna$cluster),
    clustermiRNA = as.factor(kmeans_k3_mirna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mirna' object not found. Skipping miRNA annotation.\n")
  mirna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11d: Join them all using the patient_id "glue" ---
# Start with the CNV data (the "left" side)
annotation_df_cnv <- left_join(cnv_anno_df, mrna_anno_df, by = "patient_id")
annotation_df_cnv <- left_join(annotation_df_cnv, mirna_anno_df, by = "patient_id")

# --- Step 11e: Set row names (required for pheatmap) ---
rownames(annotation_df_cnv) <- annotation_df_cnv$patient_id

cat("\n--- Merged Annotation Data Frame for Heatmaps (Head) ---\n")
print(head(annotation_df_cnv))


# ------------------------------------------------------------
# 12. Generate Heatmaps
# ------------------------------------------------------------
# We need to use the scaled, but NOT-transposed, clean matrix
# for the heatmap: 'cnv_mat_top_scaled_clean'
# And the new 'annotation_df_cnv'
# And the new 'all_kmeans_patient_ids_cnv'

# We must filter the matrix to match the (potentially smaller) patient list
# from 'data_for_kmeans_cnv'
final_patients_for_heatmap <- rownames(data_for_kmeans_cnv)
final_genes_for_heatmap <- colnames(data_for_kmeans_cnv)

heatmap_matrix <- cnv_mat_top_scaled_clean[final_genes_for_heatmap, final_patients_for_heatmap]
heatmap_annotation <- annotation_df_cnv[final_patients_for_heatmap, , drop = FALSE]

my_heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# --- Heatmap 1: K=2 ---
cat("\nGenerating k=2 heatmap (CNV)...\n")
col_order_k2 <- order(heatmap_annotation$KMeans_CNV_k2)
mat_ordered_k2 <- heatmap_matrix[, col_order_k2]
anno_ordered_k2 <- heatmap_annotation[col_order_k2, , drop = FALSE]

pheatmap(
  mat_ordered_k2,
  scale = "none", # Already scaled
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC CNV: K-Means (k=2) Clusters (Top 500 Var Genes)"
)

# --- Heatmap 2: K=3 (Ordered by Cluster 1, 2, 3) ---
cat("Generating k=3 heatmap (CNV)...\n")

# Use the sorted list of patient IDs
col_order_k3 <- all_kmeans_patient_ids_cnv 
mat_ordered_k3 <- heatmap_matrix[, col_order_k3]
anno_ordered_k3 <- heatmap_annotation[col_order_k3, , drop = FALSE]

pheatmap(
  mat_ordered_k3,
  scale = "none", # Already scaled
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k3,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC CNV: K-Means (k=3) Clusters (Top 500 Var Genes)"
)

cat("\n--- All CNV tasks complete. Check your R Plots pane. ---\n")










dataset_mut <- "mc3_gene_level/LIHC_mc3_gene_level.txt"
hub_mut <- "https://tcga.xenahubs.net" 

cat("--- Downloading Somatic Mutation (MC3) data... ---\n")

xe_mut <- XenaGenerate() %>%
  XenaFilter(filterDatasets = dataset_mut) %>%
  XenaFilter(filterHosts = hub_mut) %>%
  XenaQuery() %>%
  XenaDownload(download_R = FALSE) # Set to FALSE, we will read with fread

# ------------------------------------------------------------
# 4. Load Mutation matrix into R
# ------------------------------------------------------------
# Find the downloaded file path
filepath_mut <- list.files(
  path = tempdir(),
  pattern = "LIHC_mc3_gene_level\\.txt\\.gz$", # <-- CHANGED PATTERN
  full.names = TRUE,
  recursive = TRUE
)

cat(paste("Reading file:", filepath_mut[1], "\n"))

# Read the matrix using data.table for speed
mut_df <- fread(filepath_mut[1], sep = "\t", header = TRUE)

# ------------------------------------------------------------
# 5. Wrangle the Mutation matrix
# ------------------------------------------------------------
# This dataset has 363 samples
id_col_name <- colnames(mut_df)[1]

# 1. Check for and remove rows with NA gene IDs
na_rows <- is.na(mut_df[[id_col_name]])
if (any(na_rows)) {
  cat(paste("Found and removed", sum(na_rows), "rows with missing gene IDs.\n"))
  mut_df <- mut_df[!na_rows, ] 
}

# 2. Check for and remove rows with DUPLICATE gene IDs
dup_rows <- duplicated(mut_df[[id_col_name]])
if (any(dup_rows)) {
  cat(paste("Found and removed", sum(dup_rows), "rows with DUPLICATE gene IDs.\n"))
  mut_df <- mut_df[!dup_rows, ]
}
# --- END OF FIX ---

mut_cols <- colnames(mut_df)
mut_cols_clean <- mut_cols
mut_cols_clean[-1] <- substr(mut_cols[-1], 1, 12) # Clean column names
colnames(mut_df) <- mut_cols_clean

# Convert to a numeric matrix with gene symbols as rownames
mut_mat <- as.data.frame(mut_df)
rownames(mut_mat) <- mut_mat[[1]] # Set gene IDs as rownames
mut_mat[[1]] <- NULL             # Remove the ID column
mut_mat <- as.matrix(mut_mat)     # Convert to a numeric matrix

cat("--- Original Mutation Matrix Dimensions ---\n")
print(dim(mut_mat)) # Should be 40544 genes x 363 samples

# ------------------------------------------------------------
# 6. Filter matrix to *only* your 373 patients
# ------------------------------------------------------------
patients_to_keep_mut <- intersect(colnames(mut_mat), all_patient_ids)

cat(paste("\nFound", length(patients_to_keep_mut), "of your", 
          length(all_patient_ids), "patients in the Mutation data.\n"))

# Subset the matrix
mut_mat_filtered <- mut_mat[, patients_to_keep_mut]

cat("--- Filtered Mutation Matrix Dimensions ---\n")
print(dim(mut_mat_filtered)) # ~40544 genes x ~363 patients

cat("\n--- Head of filtered Mutation matrix ---\n")
print(mut_mat_filtered[1:4, 1:4])


# ------------------------------------------------------------
# 7. Prepare Mutation data for clustering
# ------------------------------------------------------------
cat("--- Preparing Mutation data for clustering ---\n")

# --- Step 7a: Select top 500 most variable genes ---
# We still use variance. For 0/1 data, a gene with 0.5 (50%)
# mutation rate has the highest variance.
cat("Calculating variance to find most mutated genes...\n")
mut_vars <- apply(mut_mat_filtered, 1, var, na.rm = TRUE)

# Filter for genes that *actually have variance* (var > 0)
mut_vars_clean <- mut_vars[!is.na(mut_vars) & mut_vars > 0]
mut_vars_sorted <- sort(mut_vars_clean, decreasing = TRUE)
top_500_genes_mut <- names(head(mut_vars_sorted, 500))
mut_mat_top <- mut_mat_filtered[top_500_genes_mut, ]

cat(paste("Selected top", length(top_500_genes_mut), "variable (mutated) genes.\n"))

# --- Step 7b: CRITICAL - DO NOT Z-SCORE SCALING ---
# We do NOT scale binary data.
# Instead, we just handle NAs. K-means cannot handle NAs.
# We will impute any NA to 0 (Wild Type).
if (any(is.na(mut_mat_top))) {
  cat("Found NA values... imputing to 0 (Wild Type).\n")
  mut_mat_top[is.na(mut_mat_top)] <- 0
} else {
  cat("No NA values found in top 500 genes.\n")
}

# --- Step 7c: Transpose for k-means ---
# We use the raw, unscaled 0/1 matrix
data_for_kmeans_mut <- t(mut_mat_top)

cat(paste("\nTransposed matrix for k-means. Final dimensions:",
          dim(data_for_kmeans_mut)[1], "patients x", 
          dim(data_for_kmeans_mut)[2], "genes\n"))


# ------------------------------------------------------------
# 8. Calculate Optimal 'k' (Scores)
# ------------------------------------------------------------
set.seed(123)
cat("--- Calculating Optimal 'k' (k=2 to 6) for Mutation ---\n")

# --- Elbow Method (WSS) ---
cat("Generating Elbow method plot...\n")
p_elbow_mut <- fviz_nbclust(
  data_for_kmeans_mut, 
  kmeans, 
  method = "wss",
  k.max = 6,
  nstart = 25
) + ggtitle("Mutation K-Means: Elbow Method")
print(p_elbow_mut) # Print to screen


# --- Silhouette Method ---
cat("Generating Silhouette method plot...\n")
p_sil_mut <- fviz_nbclust(
  data_for_kmeans_mut,
  kmeans,
  method = "silhouette",
  k.max = 6,
  nstart = 25
) + ggtitle("Mutation K-Means: Silhouette Method")
print(p_sil_mut) # Print to screen

cat("Optimal 'k' plots displayed.\n")

# ------------------------------------------------------------
# 9. Run K-Means for k=2 and k=3
# ------------------------------------------------------------
set.seed(123)
cat("--- Running K-Means for k=2 and k=3 ---\n")
kmeans_k2_mut <- kmeans(data_for_kmeans_mut, centers = 2, nstart = 25)
kmeans_k3_mut <- kmeans(data_for_kmeans_mut, centers = 3, nstart = 25)
cat("K-Means clustering complete.\n")

# ------------------------------------------------------------
# 10. Get Patient IDs for k=3
# ------------------------------------------------------------
cat("\n\n--- PATIENT IDs FOR K=3 CLUSTERS (Mutation) ---\n")
cluster_assignments_k3_mut <- kmeans_k3_mut$cluster
all_patient_names_mut <- names(cluster_assignments_k3_mut)

kmeans_cluster1_ids_mut <- all_patient_names_mut[cluster_assignments_k3_mut == 1]
kmeans_cluster2_ids_mut <- all_patient_names_mut[cluster_assignments_k3_mut == 2]
kmeans_cluster3_ids_mut <- all_patient_names_mut[cluster_assignments_k3_mut == 3]

all_kmeans_patient_ids_mut <- c(kmeans_cluster1_ids_mut, 
                                kmeans_cluster2_ids_mut, 
                                kmeans_cluster3_ids_mut)

# ------------------------------------------------------------
# 11. Create Merged Annotation Data Frame
# ------------------------------------------------------------
# This block merges all 4 data types (mRNA, miRNA, CNV, MUT)
# Assumes 'kmeans_k3_mrna', 'kmeans_k3_mirna', 'kmeans_k3_cnv' exist
cat("--- Building Merged Annotation Data Frame (All Data Types) ---\n")

# --- Step 11a: Create the NEW Mutation annotation data ---
mut_anno_df <- data.frame(
  patient_id = all_patient_names_mut, # Use patient_id as a COLUMN
  KMeans_MUT_k2 = as.factor(kmeans_k2_mut$cluster),
  KMeans_MUT_k3 = as.factor(kmeans_k3_mut$cluster)
)

# --- Step 11b: Create df for OLD mRNA clusters ---
if (exists("kmeans_k3_mrna")) {
  mrna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mrna$cluster),
    clustermrna = as.factor(kmeans_k3_mrna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mrna' object not found. Skipping mRNA annotation.\n")
  mrna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11c: Create df for OLD miRNA clusters ---
if (exists("kmeans_k3_mirna")) {
  mirna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mirna$cluster),
    clustermiRNA = as.factor(kmeans_k3_mirna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mirna' object not found. Skipping miRNA annotation.\n")
  mirna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11d: Create df for OLD CNV clusters ---
if (exists("kmeans_k3_cnv")) {
  cnv_anno_df <- data.frame(
    patient_id = names(kmeans_k3_cnv$cluster),
    clusterCNV = as.factor(kmeans_k3_cnv$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_cnv' object not found. Skipping CNV annotation.\n")
  cnv_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11e: Join them all using the patient_id "glue" ---
# Start with the new Mutation data (the "left" side)
annotation_df_mut <- left_join(mut_anno_df, mrna_anno_df, by = "patient_id")
annotation_df_mut <- left_join(annotation_df_mut, mirna_anno_df, by = "patient_id")
annotation_df_mut <- left_join(annotation_df_mut, cnv_anno_df, by = "patient_id")

# --- Step 11f: Set row names (required for pheatmap) ---
rownames(annotation_df_mut) <- annotation_df_mut$patient_id

cat("\n--- Merged Annotation Data Frame for Heatmaps (Head) ---\n")
print(head(annotation_df_mut))


# ------------------------------------------------------------
# 12. Generate Heatmaps
# ------------------------------------------------------------
# Use the un-scaled matrix: 'mut_mat_top'
# and the new 'annotation_df_mut'
heatmap_matrix <- mut_mat_top[colnames(data_for_kmeans_mut), rownames(data_for_kmeans_mut)]
heatmap_annotation <- annotation_df_mut[rownames(data_for_kmeans_mut), , drop = FALSE]

# --- CRITICAL CHANGE: New Color Palette ---
# Create a simple 2-color palette for 0 (WT) and 1 (Mutated)
# We set breaks to ensure 0 maps to grey and 1 maps to black
my_heatmap_colors <- c("grey95", "black")
my_heatmap_breaks <- c(-0.1, 0.5, 1.1) # 0 falls in 1st bin, 1 in 2nd

# --- Heatmap 1: K=2 ---
cat("\nGenerating k=2 heatmap (Mutation)...\n")
col_order_k2 <- order(heatmap_annotation$KMeans_MUT_k2)
mat_ordered_k2 <- heatmap_matrix[, col_order_k2]
anno_ordered_k2 <- heatmap_annotation[col_order_k2, , drop = FALSE]

pheatmap(
  mat_ordered_k2,
  scale = "none", # DATA IS NOT SCALED
  color = my_heatmap_colors,
  breaks = my_heatmap_breaks,
  annotation_col = anno_ordered_k2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC Mutation: K-Means (k=2) Clusters (Top 500 Var Genes)"
)

# --- Heatmap 2: K=3 (Ordered by Cluster 1, 2, 3) ---
cat("Generating k=3 heatmap (Mutation)...\n")

col_order_k3 <- all_kmeans_patient_ids_mut # Use the sorted list
mat_ordered_k3 <- heatmap_matrix[, col_order_k3]
anno_ordered_k3 <- heatmap_annotation[col_order_k3, , drop = FALSE]

pheatmap(
  mat_ordered_k3,
  scale = "none", # DATA IS NOT SCALED
  color = my_heatmap_colors,
  breaks = my_heatmap_breaks,
  annotation_col = anno_ordered_k3,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC Mutation: K-Means (k=3) Clusters (Top 500 Var Genes)"
)

cat("\n--- All Mutation tasks complete. Check your R Plots pane. ---\n")










# ------------------------------------------------------------
# 1. Load All Libraries (Run This Once)
# ------------------------------------------------------------
cat("--- Loading Libraries ---\n")
library(dplyr)
library(data.table)
library(stringr)
library(pheatmap)
library(factoextra)
library(ggplot2)
library(RColorBrewer)
library(XenaR)
library(UCSCXenaTools)
# ------------------------------------------------------------
# 2. Define "Ground Truth" Patient Lists
# ------------------------------------------------------------
# This assumes 'all_patient_ids' (373 total) exists from your
# first script block (the cluster_..._raw lists).
if (!exists("all_patient_ids")) {
  stop("Error: 'all_patient_ids' not found. 
       Please run the code block that defines your 3 patient clusters.")
}
cat(paste("Using master list of", length(all_patient_ids), "patients.\n"))


# ------------------------------------------------------------
# 3. Query and Download the **DNA Methylation** data
# ------------------------------------------------------------
# --- THIS IS THE MODIFIED SECTION ---

dataset_meth <- "TCGA.LIHC.sampleMap/HumanMethylation450"
hub_meth <- "https://tcga.xenahubs.net" 

cat("--- Downloading DNA Methylation data (this file is large)... ---\n")

xe_meth <- XenaGenerate() %>%
  XenaFilter(filterDatasets = dataset_meth) %>%
  XenaFilter(filterHosts = hub_meth) %>%
  XenaQuery() %>%
  XenaDownload(download_R = FALSE) # Set to FALSE, we will read with fread

# ------------------------------------------------------------
# 4. Load Methylation matrix into R
# ------------------------------------------------------------
# Find the downloaded file path
filepath_meth <- list.files(
  path = tempdir(),
  pattern = "HumanMethylation450\\.gz$", # <-- CHANGED PATTERN
  full.names = TRUE,
  recursive = TRUE
)

cat(paste("Reading file:", filepath_meth[1], "\n"))

# Read the matrix using data.table for speed
meth_df <- fread(filepath_meth[1], sep = "\t", header = TRUE)

# ------------------------------------------------------------
# 5. Wrangle the Methylation matrix (Robust Version)
# ------------------------------------------------------------
# --- START OF FIX (from Mutation script) ---
# We must clean the data frame *before* converting to a matrix
id_col_name <- colnames(meth_df)[1]

# 1. Check for and remove rows with NA probe IDs
na_rows <- is.na(meth_df[[id_col_name]])
if (any(na_rows)) {
  cat(paste("Found and removed", sum(na_rows), "rows with missing probe IDs.\n"))
  meth_df <- meth_df[!na_rows, ] 
}

# 2. Check for and remove rows with DUPLICATE probe IDs
dup_rows <- duplicated(meth_df[[id_col_name]])
if (any(dup_rows)) {
  cat(paste("Found and removed", sum(dup_rows), "rows with DUPLICATE probe IDs.\n"))
  meth_df <- meth_df[!dup_rows, ]
}
# --- END OF FIX ---

# This dataset has 429 samples
meth_cols <- colnames(meth_df)
meth_cols_clean <- meth_cols
meth_cols_clean[-1] <- substr(meth_cols[-1], 1, 12) # Clean column names
colnames(meth_df) <- meth_cols_clean

# Convert to a numeric matrix with probe symbols as rownames
meth_mat <- as.data.frame(meth_df)
rownames(meth_mat) <- meth_mat[[1]] # Set probe IDs as rownames
meth_mat[[1]] <- NULL             # Remove the ID column
meth_mat <- as.matrix(meth_mat)     # Convert to a numeric matrix

cat("--- Original (Cleaned) Methylation Matrix Dimensions ---\n")
print(dim(meth_mat)) # ~485k probes x 429 samples

# ------------------------------------------------------------
# 6. Filter matrix to *only* your 373 patients
# ------------------------------------------------------------
patients_to_keep_meth <- intersect(colnames(meth_mat), all_patient_ids)

cat(paste("\nFound", length(patients_to_keep_meth), "of your", 
          length(all_patient_ids), "patients in the Methylation data.\n"))

# Subset the matrix
meth_mat_filtered <- meth_mat[, patients_to_keep_meth]

cat("--- Filtered Methylation Matrix Dimensions ---\n")
print(dim(meth_mat_filtered)) # ~485k probes x ~373 patients

cat("\n--- Head of filtered Methylation matrix ---\n")
print(meth_mat_filtered[1:4, 1:4])


# ------------------------------------------------------------
# 7. Prepare Methylation data for clustering (Robust Version)
# ------------------------------------------------------------
cat("--- Preparing Methylation data for clustering ---\n")

# --- Step 7a: Select top 500 most variable probes ---
cat("Calculating variance (this may take a moment)...\n")
meth_vars <- apply(meth_mat_filtered, 1, var, na.rm = TRUE)
meth_vars_clean <- meth_vars[!is.na(meth_vars) & meth_vars > 0]
meth_vars_sorted <- sort(meth_vars_clean, decreasing = TRUE)
top_500_probes <- names(head(meth_vars_sorted, 500))
meth_mat_top <- meth_mat_filtered[top_500_probes, ]

cat(paste("Selected top", length(top_500_probes), "variable probes.\n"))

# --- Step 7b: Z-score the data (Scale by probe/row) ---
cat("Scaling data (Z-score by row) and ignoring NAs...\n")
row_means <- rowMeans(meth_mat_top, na.rm = TRUE)
row_sds <- apply(meth_mat_top, 1, sd, na.rm = TRUE)
meth_mat_top_scaled <- (meth_mat_top - row_means) / row_sds

# --- Step 7c: Remove problematic probes (NA/NaN/Inf) ---
# !is.finite() checks for NA, NaN, Inf
rows_with_nonfinite <- apply(meth_mat_top_scaled, 1, function(row) any(!is.finite(row)))

if(any(rows_with_nonfinite)) {
  cat(paste("Found and removed", sum(rows_with_nonfinite), "probe(s) with NA, NaN, or Inf values (likely 0 variance).\n"))
  meth_mat_top_scaled_clean <- meth_mat_top_scaled[!rows_with_nonfinite, ]
} else {
  cat("No problematic probes (NA/NaN/Inf) found.\n")
  meth_mat_top_scaled_clean <- meth_mat_top_scaled
}

# --- Step 7d: Transpose for k-means ---
data_for_kmeans_meth_raw <- t(meth_mat_top_scaled_clean)

# --- Step 7e (FINAL CHECK): Remove problematic patients ---
patients_with_nonfinite <- apply(data_for_kmeans_meth_raw, 1, function(row) any(!is.finite(row)))

if(any(patients_with_nonfinite)) {
  cat(paste("Found and removed", sum(patients_with_nonfinite), "patient(s) with NA/NaN/Inf values.\n"))
  data_for_kmeans_meth <- data_for_kmeans_meth_raw[!patients_with_nonfinite, ]
} else {
  cat("All patients have finite data. Ready for clustering.\n")
  data_for_kmeans_meth <- data_for_kmeans_meth_raw
}

cat(paste("\nTransposed matrix for k-means. Final dimensions:",
          dim(data_for_kmeans_meth)[1], "patients x", 
          dim(data_for_kmeans_meth)[2], "probes\n"))


# ------------------------------------------------------------
# 8. Calculate Optimal 'k' (Scores)
# ------------------------------------------------------------
set.seed(123)
cat("--- Calculating Optimal 'k' (k=2 to 6) for Methylation ---\n")

# --- Elbow Method (WSS) ---
cat("Generating Elbow method plot...\n")
p_elbow_meth <- fviz_nbclust(
  data_for_kmeans_meth, 
  kmeans, 
  method = "wss",
  k.max = 6,
  nstart = 25
) + ggtitle("DNA Methylation K-Means: Elbow Method")
print(p_elbow_meth) # Print to screen


# --- Silhouette Method ---
cat("Generating Silhouette method plot...\n")
p_sil_meth <- fviz_nbclust(
  data_for_kmeans_meth,
  kmeans,
  method = "silhouette",
  k.max = 6,
  nstart = 25
) + ggtitle("DNA Methylation K-Means: Silhouette Method")
print(p_sil_meth) # Print to screen

cat("Optimal 'k' plots displayed.\n")

# ------------------------------------------------------------
# 9. Run K-Means for k=2 and k=3
# ------------------------------------------------------------
set.seed(123)
cat("--- Running K-Means for k=2 and k=3 ---\n")
kmeans_k2_meth <- kmeans(data_for_kmeans_meth, centers = 2, nstart = 25)
kmeans_k3_meth <- kmeans(data_for_kmeans_meth, centers = 3, nstart = 25)
cat("K-Means clustering complete.\n")

# ------------------------------------------------------------
# 10. Get Patient IDs for k=3
# ------------------------------------------------------------
cat("\n\n--- PATIENT IDs FOR K=3 CLUSTERS (Methylation) ---\n")
cluster_assignments_k3_meth <- kmeans_k3_meth$cluster
all_patient_names_meth <- names(cluster_assignments_k3_meth)

kmeans_cluster1_ids_meth <- all_patient_names_meth[cluster_assignments_k3_meth == 1]
kmeans_cluster2_ids_meth <- all_patient_names_meth[cluster_assignments_k3_meth == 2]
kmeans_cluster3_ids_meth <- all_patient_names_meth[cluster_assignments_k3_meth == 3]

all_kmeans_patient_ids_meth <- c(kmeans_cluster1_ids_meth, 
                                 kmeans_cluster2_ids_meth, 
                                 kmeans_cluster3_ids_meth)

# ------------------------------------------------------------
# 11. Create Merged Annotation Data Frame
# ------------------------------------------------------------
# This block merges all 5 data types (mRNA, miRNA, CNV, MUT, METH)
# Assumes previous 'kmeans_k3_...' objects exist
cat("--- Building Merged Annotation Data Frame (All 5 Data Types) ---\n")

# --- Step 11a: Create the NEW Methylation annotation data ---
meth_anno_df <- data.frame(
  patient_id = all_patient_names_meth, # Use patient_id as a COLUMN
  KMeans_METH_k2 = as.factor(kmeans_k2_meth$cluster),
  KMeans_METH_k3 = as.factor(kmeans_k3_meth$cluster),
  clusterMeth = as.factor(kmeans_k3_meth$cluster) # Your requested name
)

# --- Step 11b: Create df for OLD mRNA clusters ---
if (exists("kmeans_k3_mrna")) {
  mrna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mrna$cluster),
    clustermrna = as.factor(kmeans_k3_mrna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mrna' object not found. Skipping mRNA annotation.\n")
  mrna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11c: Create df for OLD miRNA clusters ---
if (exists("kmeans_k3_mirna")) {
  mirna_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mirna$cluster),
    clustermiRNA = as.factor(kmeans_k3_mirna$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mirna' object not found. Skipping miRNA annotation.\n")
  mirna_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11d: Create df for OLD CNV clusters ---
if (exists("kmeans_k3_cnv")) {
  cnv_anno_df <- data.frame(
    patient_id = names(kmeans_k3_cnv$cluster),
    clusterCNV = as.factor(kmeans_k3_cnv$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_cnv' object not found. Skipping CNV annotation.\n")
  cnv_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11e: Create df for OLD Mutation clusters ---
if (exists("kmeans_k3_mut")) {
  mut_anno_df <- data.frame(
    patient_id = names(kmeans_k3_mut$cluster),
    clusterMut = as.factor(kmeans_k3_mut$cluster)
  )
} else {
  cat("Warning: 'kmeans_k3_mut' object not found. Skipping Mutation annotation.\n")
  mut_anno_df <- data.frame(patient_id = character(0)) # Empty df
}

# --- Step 11f: Join them all using the patient_id "glue" ---
# Start with the new Methylation data (the "left" side)
annotation_df_meth <- left_join(meth_anno_df, mrna_anno_df, by = "patient_id")
annotation_df_meth <- left_join(annotation_df_meth, mirna_anno_df, by = "patient_id")
annotation_df_meth <- left_join(annotation_df_meth, cnv_anno_df, by = "patient_id")
annotation_df_meth <- left_join(annotation_df_meth, mut_anno_df, by = "patient_id")

# --- Step 11g: Set row names (required for pheatmap) ---
rownames(annotation_df_meth) <- annotation_df_meth$patient_id

cat("\n--- Merged Annotation Data Frame for Heatmaps (Head) ---\n")
print(head(annotation_df_meth))


# ------------------------------------------------------------
# 12. Generate Heatmaps
# ------------------------------------------------------------
# Use the Z-scored matrix: 'meth_mat_top_scaled_clean'
# Filtered to the final list of patients and probes
final_patients_for_heatmap <- rownames(data_for_kmeans_meth)
final_probes_for_heatmap <- colnames(data_for_kmeans_meth)

heatmap_matrix <- meth_mat_top_scaled_clean[final_probes_for_heatmap, final_patients_for_heatmap]
heatmap_annotation <- annotation_df_meth[final_patients_for_heatmap, , drop = FALSE]

# --- Color Palette: Back to Blue-Red for Z-scored data ---
my_heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# --- Heatmap 1: K=2 ---
cat("\nGenerating k=2 heatmap (Methylation)...\n")
col_order_k2 <- order(heatmap_annotation$KMeans_METH_k2)
mat_ordered_k2 <- heatmap_matrix[, col_order_k2]
anno_ordered_k2 <- heatmap_annotation[col_order_k2, , drop = FALSE]

pheatmap(
  mat_ordered_k2,
  scale = "none", # Already scaled
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC Methylation: K-Means (k=2) Clusters (Top 500 Var Probes)"
)

# --- Heatmap 2: K=3 (Ordered by Cluster 1, 2, 3) ---
cat("Generating k=3 heatmap (Methylation)...\n")

col_order_k3 <- all_kmeans_patient_ids_meth # Use the sorted list
mat_ordered_k3 <- heatmap_matrix[, col_order_k3]
anno_ordered_k3 <- heatmap_annotation[col_order_k3, , drop = FALSE]

pheatmap(
  mat_ordered_k3,
  scale = "none", # Already scaled
  color = my_heatmap_colors,
  annotation_col = anno_ordered_k3,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "TCGA-LIHC Methylation: K-Means (k=3) Clusters (Top 500 Var Probes)"
)

cat("\n--- All Methylation tasks complete. Check your R Plots pane. ---\n")






# ------------------------------------------------------------
# 1. Install and Load Libraries
# ------------------------------------------------------------
cat("--- Installing/Loading Libraries ---\n")

# Install MOVICS from CRAN
if (!require(MOVICS)) install.packages("MOVICS")

# MOVICS has many dependencies. If the above fails, you may need
# BiocManager to install specific packages, e.g.:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("ConsensusClusterPlus", "iClusterPlus", "survival"))

# Load all required libraries
library(MOVICS)
library(dplyr) # For data manipulation

# ------------------------------------------------------------
# 2. CRITICAL: Prepare Data for Multi-Omics Integration
# ------------------------------------------------------------
# MOVICS requires all data matrices to have:
# 1. Patients in ROWS
# 2. Features in COLUMNS
# 3. The *exact same set of patients* in the *exact same order*.

# Your 'data_for_kmeans_...' objects already have patients in rows,
# but they have different patient lists due to separate NA cleaning.
# We must find the common patients *across all 5 datasets*.

cat("--- Finding Common Patients Across All 5 Datasets ---\n")

# Get the patient IDs (rownames) from each dataset
patients_mrna <- rownames(data_for_kmeans_mrna)
patients_mirna <- rownames(data_for_kmeans_mirna)
patients_meth <- rownames(data_for_kmeans_meth)
patients_cnv <- rownames(data_for_kmeans_cnv)
patients_mut <- rownames(data_for_kmeans_mut)

# Create a list of all patient vectors
patient_lists <- list(
  mrna = patients_mrna,
  mirna = patients_mirna,
  meth = patients_meth,
  cnv = patients_cnv,
  mut = patients_mut
)

# Find the intersection (common patients)
common_patients <- Reduce(intersect, patient_lists)

cat(paste("Found", length(common_patients), "common patients across all 5 datasets.\n"))
print(head(data_for_kmeans_mrna))
# --- Create the final list of data matrices ---
# Filter each matrix to only these common patients, in the same order.
mo.data.list <- list(
  mRNA = data_for_kmeans_mrna[common_patients, ],
  miRNA = data_for_kmeans_mirna[common_patients, ],
  Methylation = data_for_kmeans_meth[common_patients, ],
  CNV = data_for_kmeans_cnv[common_patients, ],
  Mutation = data_for_kmeans_mut[common_patients, ]
)

omics_to_transpose <- mo.data.list[c("mRNA", "miRNA", "Methylation", "CNV")]

transposed_omics <- lapply(omics_to_transpose, t)
print(sapply(mo.data.list["Mutation"], dim))

# 2. Combine the transposed matrices with the Mutation matrix (which is not transposed)
#    We use append() and include the Mutation matrix directly.
mo.data.final <- append(transposed_omics, 
                        lapply(mo.data.list["Mutation"], t))


# 3. Reorder the final list back to the original sequence (mRNA...Mutation)
mo.data.final <- mo.data.final[names(mo.data.list)]

# --- Verification ---
cat("Final Structure (Rows x Columns):\n")
print(sapply(mo.data.final, dim))


print(head(mo.data.final))




mo.data.for.clustering <- mo.data.final

print(sapply(mo.data.for.clustering, dim))

is.binary.vector <- c(F, F, F, F, T) 

cat("--- Calculating optimal cluster number (getClustNum) ---\n")
optk <- getClustNum(data = mo.data.for.clustering,
                    is.binary = is.binary.vector,
                    try.N.clust = 2:8, 
                    fig.name = "OPTIMAL_CLUSTER_NUMBER")

# -------------------------------------------------------------------
#   *** 🛑 ACTION REQUIRED 🛑 ***
#
#   1. Inspect "OPTIMAL_CLUSTER_NUMBER.pdf" and set N.CLUST.
#   2. The script must be stopped here.
#
# -------------------------------------------------------------------

# You must set N.CLUST based on your plot before running the next block.
N.CLUST <- 2  # <<< Using k=3 as an example for demonstration >>>

# ------------------------------------------------------------
# 2. RUN ALL 10 MOIC ALGORITHMS & CONSENSUS
# ------------------------------------------------------------
cat(paste("--- Running 10 MOIC algorithms for k =", N.CLUST, "... ---\n"))

# Define type vector for 5 omics: (G, G, G, G, B)
data.types.vector <- c("gaussian", "gaussian", "gaussian", "gaussian", "binomial")

moic.res.list <- getMOIC(
  data = mo.data.for.clustering,
  methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster",
                     "ConsensusClustering", "IntNMF", "CIMLR",
                     "MoCluster", "iClusterBayes"),
  N.clust = N.CLUST,
  type = data.types.vector
)

# Combine results into one robust clustering
cmoic.res <- getConsensusMOIC(
  moic.res.list = moic.res.list,
  fig.name = "CONSENSUS_HEATMAP",
  distance = "euclidean",
  linkage = "average"
)

cat("--- Consensus Clustering Complete. ---\n")


# ------------------------------------------------------------
# 3. STANDARDIZATION & PLOTDATA PREP
# ------------------------------------------------------------
cat("--- Standardizing Data and M-Value Conversion ---\n")

# Copy the clustered data (Features x Patients)
indata <- mo.data.for.clustering 


# 3.2 Standardization (Z-score and Truncation)
# Vectors adapted for 5 omics (F, F, F, F, T)
plotdata <- getStdiz(data = indata,
                     halfwidth = c(2, 2, 2, 2, NA), # NA for Mutation
                     centerFlag = c(T, T, T, T, F), 
                     scaleFlag = c(T, T, T, T, F))  

cat("--- Plotting Data Standardization Complete. ---\n")


# ------------------------------------------------------------
# 4. FINAL PLOTTING: Comprehensive Heatmap
# ------------------------------------------------------------
cat("--- Generating FINAL COMPREHENSIVE HEATMAP ---\n")

# --- Define Annotation Objects (using placeholders) ---
# NOTE: Replace with your final annotation objects if needed.
dataset_names <- names(mo.data.for.clustering)
feat <- moic.res.list$iClusterBayes$feat.res # Use iClusterBayes for feature ranking
annRow <- list(
  feat[which(feat$dataset == dataset_names[1]),][1:10,"feature"], # mRNA
  feat[which(feat$dataset == dataset_names[2]),][1:10,"feature"], # miRNA
  feat[which(feat$dataset == dataset_names[3]),][1:10,"feature"], # Methylation
  feat[which(feat$dataset == dataset_names[4]),][1:10,"feature"], # CNV
  feat[which(feat$dataset == dataset_names[5]),][1:10,"feature"]  # Mutation
)
names(annRow) <- dataset_names

# Define 5-omics color list
col.list <- list(
  c("#00FF00", "black", "#FF0000"), # mRNA (Green-Black-Red)
  c("#6699CC", "white", "#FF3C38"), # miRNA
  c("#0074FE", "white", "#F00003"), # Methylation (Blue-White-Red)
  c("darkblue", "white", "darkred"), # CNV
  c("grey90", "black") # Mutation (Binary)
)


getMoHeatmap(data = plotdata, # Use the standardized plotting data (Features x Patients)
             row.title = c("mRNA Expression","miRNA Expression","DNA Methylation","CNV Status","Mutation Status"),
             is.binary = c(F, F, F, F, T),
             legend.name = c("mRNA Z-score","miRNA Z-score","M-value Z-score","CNV Z-score","Mutated"),
             
             clust.res = cmoic.res$clust.res, # Consensus result
             clust.dend = cmoic.res$clust.dend,
             annRow = annRow, # Annotate with top features
             color = col.list, 
             
             # General Plot Settings
             show.rownames = c(F, F, F, F, F), 
             show.colnames = FALSE,
             annCol = NULL, 
             annColors = NULL, 
             fig.name = "COMPREHENSIVE_HEATMAP_CONSENSUS")

cat("\n✅ Full Clustering and Heatmap Generation Complete! You can now proceed to the COMP and RUN Modules for characterization.")

cat("\n✅ Full Clustering and Heatmap Generation Complete! You can now proceed to the COMP and RUN Modules for characterization.")



cat("\n=================================================================\n")
cat(paste("=== PATIENT SUBTYPES BASED ON CONSENSUS CLUSTERING (k =", N.CLUST, ") ===\n"))
cat("=================================================================\n")

# Extract the final clustering result data frame
cluster_assignments <- cmoic.res$clust.res

# Ensure the cluster column is treated as a factor for correct grouping
cluster_assignments$clust <- as.factor(cluster_assignments$clust)

# Sort patients by their cluster assignment
cluster_assignments_sorted <- cluster_assignments[order(cluster_assignments$clust), ]

# Loop through each unique cluster and print the patient IDs
for (k in sort(unique(cluster_assignments_sorted$clust))) {
  
  # Filter patients belonging to the current cluster
  patients_in_cluster <- cluster_assignments_sorted$samID[cluster_assignments_sorted$clust == k]
  
  cluster_label <- paste0("CS", k) # Create a standardized label (e.g., CS1, CS2, CS3)
  
  cat("\n----------------------------------")
  cat(paste("\nSubtype:", cluster_label))
  cat(paste("\nTotal Patients:", length(patients_in_cluster)))
  cat("\n----------------------------------\n")
  
  # CRITICAL FIX: Use cat() with sep="\n" to print each ID on a new line 
  # without vector brackets or quotes.
  cat(patients_in_cluster, sep = "\n")
  
}


# You must set N.CLUST based on your plot before running the next block.
N.CLUST <- 3  # <<< Using k=3 as an example for demonstration >>>

# ------------------------------------------------------------
# 2. RUN ALL 10 MOIC ALGORITHMS & CONSENSUS
# ------------------------------------------------------------
cat(paste("--- Running 10 MOIC algorithms for k =", N.CLUST, "... ---\n"))

# Define type vector for 5 omics: (G, G, G, G, B)
data.types.vector <- c("gaussian", "gaussian", "gaussian", "gaussian", "binomial")

moic.res.list <- getMOIC(
  data = mo.data.for.clustering,
  methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster",
                     "ConsensusClustering", "IntNMF", "CIMLR",
                     "MoCluster", "iClusterBayes"),
  N.clust = N.CLUST,
  type = data.types.vector
)

# Combine results into one robust clustering
cmoic.res <- getConsensusMOIC(
  moic.res.list = moic.res.list,
  fig.name = "CONSENSUS_HEATMAP",
  distance = "euclidean",
  linkage = "average"
)

cat("--- Consensus Clustering Complete. ---\n")


# ------------------------------------------------------------
# 3. STANDARDIZATION & PLOTDATA PREP
# ------------------------------------------------------------
cat("--- Standardizing Data and M-Value Conversion ---\n")

# Copy the clustered data (Features x Patients)
indata <- mo.data.for.clustering 


# 3.2 Standardization (Z-score and Truncation)
# Vectors adapted for 5 omics (F, F, F, F, T)
plotdata <- getStdiz(data = indata,
                     halfwidth = c(2, 2, 2, 2, NA), # NA for Mutation
                     centerFlag = c(T, T, T, T, F), 
                     scaleFlag = c(T, T, T, T, F))  

cat("--- Plotting Data Standardization Complete. ---\n")


# ------------------------------------------------------------
# 4. FINAL PLOTTING: Comprehensive Heatmap
# ------------------------------------------------------------
cat("--- Generating FINAL COMPREHENSIVE HEATMAP ---\n")

# --- Define Annotation Objects (using placeholders) ---
# NOTE: Replace with your final annotation objects if needed.
dataset_names <- names(mo.data.for.clustering)
feat <- moic.res.list$iClusterBayes$feat.res # Use iClusterBayes for feature ranking
annRow <- list(
  feat[which(feat$dataset == dataset_names[1]),][1:10,"feature"], # mRNA
  feat[which(feat$dataset == dataset_names[2]),][1:10,"feature"], # miRNA
  feat[which(feat$dataset == dataset_names[3]),][1:10,"feature"], # Methylation
  feat[which(feat$dataset == dataset_names[4]),][1:10,"feature"], # CNV
  feat[which(feat$dataset == dataset_names[5]),][1:10,"feature"]  # Mutation
)
names(annRow) <- dataset_names

# Define 5-omics color list
col.list <- list(
  c("#00FF00", "black", "#FF0000"), # mRNA (Green-Black-Red)
  c("#6699CC", "white", "#FF3C38"), # miRNA
  c("#0074FE", "white", "#F00003"), # Methylation (Blue-White-Red)
  c("darkblue", "white", "darkred"), # CNV
  c("grey90", "black") # Mutation (Binary)
)


getMoHeatmap(data = plotdata, # Use the standardized plotting data (Features x Patients)
             row.title = c("mRNA Expression","miRNA Expression","DNA Methylation","CNV Status","Mutation Status"),
             is.binary = c(F, F, F, F, T),
             legend.name = c("mRNA Z-score","miRNA Z-score","M-value Z-score","CNV Z-score","Mutated"),
             
             clust.res = cmoic.res$clust.res, # Consensus result
             clust.dend = cmoic.res$clust.dend,
             annRow = annRow, # Annotate with top features
             color = col.list, 
             
             # General Plot Settings
             show.rownames = c(F, F, F, F, F), 
             show.colnames = FALSE,
             annCol = NULL, 
             annColors = NULL, 
             fig.name = "COMPREHENSIVE_HEATMAP_CONSENSUS")

cat("\n✅ Full Clustering and Heatmap Generation Complete! You can now proceed to the COMP and RUN Modules for characterization.")



cat("\n=================================================================\n")
cat(paste("=== PATIENT SUBTYPES BASED ON CONSENSUS CLUSTERING (k =", N.CLUST, ") ===\n"))
cat("=================================================================\n")

# Extract the final clustering result data frame
cluster_assignments <- cmoic.res$clust.res

# Ensure the cluster column is treated as a factor for correct grouping
cluster_assignments$clust <- as.factor(cluster_assignments$clust)

# Sort patients by their cluster assignment
cluster_assignments_sorted <- cluster_assignments[order(cluster_assignments$clust), ]

