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
names(brca.tcga)
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
  color = my_blue_red_colors,
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

#-------------------


# --- 0. Prerequisite ---
# Assumes 'all_patient_ids' (your list of 373 patients) is in the environment.
# Load libraries
library(pheatmap)
library(clValid)
library(data.table)

# --- 1. Download Copy Number Data ---
cat("Downloading TCGA-LIHC.gene-level_ascat2.tsv.gz...\n")
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.gene-level_ascat2.tsv.gz"
dest_file <- "TCGA-LIHC.gene-level_ascat2.tsv.gz"

if (!file.exists(dest_file)) {
  download.file(url, dest_file)
}

# --- 2. Load Data ---
cnv_data <- read.delim(gzfile(dest_file), 
                       row.names = 1, 
                       check.names = FALSE)
# Clean column names to 12-char ID
colnames(cnv_data) <- substr(colnames(cnv_data), 1, 12)
cat("Download complete.\n")

# --- 3. Filter for your Patient Cohort (NEW STEP) ---
# Find the patients from your list that are in this dataset
common_patients_cna <- intersect(colnames(cnv_data), all_patient_ids)
cat(paste("Found", length(common_patients_cna), "of your patients in the CNA data.\n"))

# Filter the matrix to ONLY these patients
cnv_data_filtered <- cnv_data[, common_patients_cna]

# --- 4. Data Preparation (on Filtered Data) ---

### Find Top 500 Variable Genes
# Calculate variance for each gene (row) *using the filtered data*
row_vars <- apply(cnv_data_filtered, 1, var)
row_vars <- row_vars[!is.na(row_vars) & row_vars > 0]
top_500_genes <- names(sort(row_vars, decreasing = TRUE)[1:500])
data_subset <- cnv_data_filtered[top_500_genes, ]

### Row-Normalize Data (Z-score)
scaled_data <- t(scale(t(data_subset)))
scaled_data[is.na(scaled_data)] <- 0

# --- 5. Heatmap Visualization (k=2 and k=3) ---

cat("Generating heatmap for k=3...\n")
pheatmap(
  scaled_data,
  main = "Top 500 Variable Genes (CNA, k=3, Filtered Cohort)",
  scale = "none",
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  cutree_cols = 3 
)

# --- 6. Internal Cluster Validation (k=2 to 6) ---

cat("Running internal cluster validation for k=2 to 6...\n")
cl_validation_cna <- clValid(
  obj = t(scaled_data), 
  nClust = 2:6,
  clMethods = "kmeans",
  validation = "internal"
)

cat("--- Cluster Validation Results (CNA) ---\n")
print(summary(cl_validation_cna))

# --- 1. Define the Data for Clustering ---
# This is the scaled, top 500 gene matrix from your CNV script
# We transpose it so patients are rows
data_for_kmeans_cnv <- t(scaled_data)

# --- 1. Define the Data for Clustering ---
# This is the scaled, top 500 gene matrix from your CNV script
# We transpose it so patients are rows
data_for_kmeans_cnv <- t(scaled_data)

# Set seed for reproducible results
set.seed(123) 

# --- 2. Run k-means for k=2 AND k=3 ---
kmeans_k2_cnv <- kmeans(data_for_kmeans_cnv, centers = 2, nstart = 25)
kmeans_k3_cnv <- kmeans(data_for_kmeans_cnv, centers = 3, nstart = 25)

# Set seed for reproducible results
# --- 3. Create Comprehensive Annotation DF ---
# Get all patient IDs that are in your final CNV matrix
all_patients_in_cnv_matrix <- colnames(scaled_data)

# Create the data frame
cnv_annotation_df <- data.frame(
  row.names = all_patients_in_cnv_matrix
)

# --- 3a. Add Original 'Cluster' assignments ---
# (Assumes cluster1_ids, cluster2_ids, cluster3_ids exist)
c1_in_cnv <- intersect(all_patients_in_cnv_matrix, cluster1_ids)
c2_in_cnv <- intersect(all_patients_in_cnv_matrix, cluster2_ids)
c3_in_cnv <- intersect(all_patients_in_cnv_matrix, cluster3_ids)

cnv_annotation_df[c1_in_cnv, "Cluster"] <- "Cluster 1"
cnv_annotation_df[c2_in_cnv, "Cluster"] <- "Cluster 2"
cnv_annotation_df[c3_in_cnv, "Cluster"] <- "Cluster 3"

# --- 3b. Add New k-means assignments ---
cnv_annotation_df$KMeans_k2 <- as.factor(kmeans_k2_cnv$cluster)
cnv_annotation_df$KMeans_k3 <- as.factor(kmeans_k3_cnv$cluster)

# Convert original cluster to a factor for correct ordering
cnv_annotation_df$Cluster <- factor(cnv_annotation_df$Cluster, 
                                    levels = c("Cluster 1", "Cluster 2", "Cluster 3"))

cat("\n--- CNV Annotation DF Head ---\n")
print(head(cnv_annotation_df))
# --- 4. Create the Nested Sorting Order ---
col_order_nested_cnv <- order(
  cnv_annotation_df$Cluster,
  cnv_annotation_df$KMeans_k3,
  cnv_annotation_df$KMeans_k2
)

# --- 5. Re-order the matrix and annotation ---
cnv_mat_top_ordered <- scaled_data[, col_order_nested_cnv]
cnv_annotation_df_ordered <- cnv_annotation_df[col_order_nested_cnv, , drop = FALSE]

# --- 6. Get Gap Locations (based on primary sort key) ---
gaps_cnv <- cumsum(table(cnv_annotation_df_ordered$Cluster))

# --- 7. Plot the Heatmap ---
# (Using the diverging color palette for Z-score data)
my_blue_red_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
pheatmap(
  cnv_mat_top_ordered,         # <-- Use the RE-ORDERED CNV matrix
  scale = "none",                # <-- Data is already Z-scored
  color = my_blue_red_colors,         # <-- Use the diverging palette
  annotation_col = cnv_annotation_df_ordered, # <-- Use the RE-ORDERED annotation
  cluster_cols = FALSE,          # <-- MUST BE FALSE to keep your order
  cluster_rows = TRUE,           # <-- Keep clustering genes
  gaps_col = gaps_cnv[-length(gaps_cnv)], # <-- Adds gaps between original clusters
  fontsize_row = 6,
  fontsize_col = 6,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "TCGA-LIHC CNV: Top 500 Genes (Grouped by Cluster, K3, K2)"
)


library(factoextra)

# --- 3. Compute Scores: Elbow Method (WSS) ---
cat("\n--- Generating Elbow Method Plot (CNV) ---\n")
fviz_nbclust(
  data_for_kmeans_cnv, 
  kmeans, 
  method = "wss", 
  k.max = 6, 
  nstart = 25 
)

# --- 4. Compute Scores: Silhouette Method ---
cat("\n--- Generating Silhouette Method Plot (CNV) ---\n")
fviz_nbclust(
  data_for_kmeans_cnv, 
  kmeans, 
  method = "silhouette",
  k.max = 6, 
  nstart = 25
)


# --- 0. Load Libraries ---
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(factoextra)

# --- 1. Load Data ---
cat("Reading data file with data.table::fread()...\n")
dest_file <- "TCGA-LIHC.methylation450.tsv.gz"
# (Add download.file() code here if you don't have the file)

meth_data_df <- fread(dest_file, 
                      data.table = FALSE, 
                      check.names = FALSE)
rownames(meth_data_df) <- meth_data_df[, 1]
meth_data <- meth_data_df[, -1]
# Clean column names to 12-char ID
colnames(meth_data) <- substr(colnames(meth_data), 1, 12)
cat("Read complete.\n")
# --- 2. Filter for your Patient Cohort ---
common_patients_meth <- intersect(colnames(meth_data), all_patient_ids)
cat(paste("Found", length(common_patients_meth), "of your patients in the Methylation data.\n"))
meth_data_filtered <- meth_data[, common_patients_meth]

# --- 3. Data Preparation ---
cat("Calculating variance for 486,000+ probes...\n")
row_vars <- apply(meth_data_filtered, 1, var, na.rm = TRUE) # Use na.rm=TRUE
row_vars <- row_vars[!is.na(row_vars) & row_vars > 0]
top_500_probes <- names(sort(row_vars, decreasing = TRUE)[1:500])
data_subset <- meth_data_filtered[top_500_probes, ]
cat("Top 500 probes selected.\n")

# Row-Normalize Data (Z-score)
scaled_data <- t(scale(t(data_subset)))
scaled_data[is.na(scaled_data)] <- 0

# --- 4. Run K-Means (k=2 & k=3) ---
data_for_kmeans_meth <- t(scaled_data) # Patients as rows
set.seed(123) 
kmeans_k2_meth <- kmeans(data_for_kmeans_meth, centers = 2, nstart = 25)
kmeans_k3_meth <- kmeans(data_for_kmeans_meth, centers = 3, nstart = 25)

# --- 5. Create Comprehensive Annotation ---
all_patients_in_meth_matrix <- colnames(scaled_data)
meth_annotation_df <- data.frame(row.names = all_patients_in_meth_matrix)

# Add Original 'Cluster'
c1_in_meth <- intersect(all_patients_in_meth_matrix, cluster1_ids)
c2_in_meth <- intersect(all_patients_in_meth_matrix, cluster2_ids)
c3_in_meth <- intersect(all_patients_in_meth_matrix, cluster3_ids)
meth_annotation_df[c1_in_meth, "Cluster"] <- "Cluster 1"
meth_annotation_df[c2_in_meth, "Cluster"] <- "Cluster 2"
meth_annotation_df[c3_in_meth, "Cluster"] <- "Cluster 3"

# Add New k-means
meth_annotation_df$KMeans_k2 <- as.factor(kmeans_k2_meth$cluster)
meth_annotation_df$KMeans_k3 <- as.factor(kmeans_k3_meth$cluster)
meth_annotation_df$Cluster <- factor(meth_annotation_df$Cluster, 
                                     levels = c("Cluster 1", "Cluster 2", "Cluster 3"))
cat("\n--- Methylation Annotation DF Head ---\n")
print(head(meth_annotation_df))

# --- 6. Generate the Sorted Heatmap ---
# Create the nested sorting order
col_order_nested_meth <- order(
  meth_annotation_df$Cluster,
  meth_annotation_df$KMeans_k3,
  meth_annotation_df$KMeans_k2
)
# Re-order the matrix and annotation
meth_mat_top_ordered <- scaled_data[, col_order_nested_meth]
meth_annotation_df_ordered <- meth_annotation_df[col_order_nested_meth, , drop = FALSE]

# Get gap locations
gaps_meth <- cumsum(table(meth_annotation_df_ordered$Cluster))

# Create the Blue-White-Red color palette
my_blue_red_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Plot
pheatmap(
  meth_mat_top_ordered,
  scale = "none",
  color = my_blue_red_colors,        # <-- Use Blue-Red palette
  annotation_col = meth_annotation_df_ordered,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  gaps_col = gaps_meth[-length(gaps_meth)],
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "TCGA-LIHC Methylation: Top 500 Probes (Grouped by Cluster, K3, K2)"
)


# --- 7. Compute Optimal Cluster Scores (factoextra) ---
cat("\n--- Generating Elbow Method Plot (Methylation) ---\n")
fviz_nbclust(
  data_for_kmeans_meth, 
  kmeans, 
  method = "wss", 
  k.max = 6, 
  nstart = 25 
)

cat("\n--- Generating Silhouette Method Plot (Methylation) ---\n")
fviz_nbclust(
  data_for_kmeans_meth, 
  kmeans, 
  method = "silhouette",
  k.max = 6, 
  nstart = 25
)


