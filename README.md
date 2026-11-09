# pangenome_SynComs
pangenome analysis of Atsphere SynComs
# PRJNA297942 Assembly Details

#beta-proteobacteria
#Step 1: Create beta_accessions_30.txt

cat << 'EOF' > beta_accessions_30.txt
GCA_001428165.1
GCA_001427425.1
GCA_001428845.1
GCA_001428645.1
GCA_001428665.1
GCA_001429215.1
GCA_001428965.1
GCA_001425125.1
GCA_001427445.1
GCA_001428225.1
GCA_001427305.1
GCA_001428585.1
GCA_001425005.1
GCA_001428625.1
GCA_001426525.1
GCA_001427325.1
GCA_001425685.1
GCA_001425045.1
GCA_001425265.1
GCA_001425905.1
GCA_001428545.1
GCA_001425705.1
GCA_001425055.1
GCA_001426885.1
GCA_001425165.1
GCA_001427705.1
GCA_001424835.1
GCA_001425205.1
GCA_001426595.1
GCA_001426505.1
EOF

#Step 2: Download and Extract the FASTA Files
# 1. Ensure Conda environment is active
conda activate ppanggolin_env

# 2. Download assemblies
datasets download genome accession --inputfile beta_accessions_30.txt \
  --include genome \
  --filename beta_genomes_30.zip

# 3. Extract and copy to fasta_input/
unzip -n beta_genomes_30.zip -d ncbi_beta_extracted
mkdir -p fasta_input
cp -n ncbi_beta_extracted/ncbi_dataset/data/*/*.fna fasta_input/

nano make_beta_list.sh
#!/bin/bash

# Initialize clean list
> beta_genomes_30.fasta.list

# Match accessions to absolute paths
while IFS= read -r acc; do
    file=$(ls fasta_input/${acc}*.fna 2>/dev/null | head -n 1)
    if [ -f "$file" ]; then
        abs_path=$(readlink -f "$file")
        echo -e "${acc}\t${abs_path}" >> beta_genomes_30.fasta.list
    else
        echo "Warning: FASTA file for ${acc} not found!"
    fi
done < beta_accessions_30.txt

echo "Done! Total entries in list:"
wc -l beta_genomes_30.fasta.list

bash make_beta_list.sh

#Step 4: Run PPanGGOLiN for the beta Group
ppanggolin all \
  --fasta beta_genomes_30.fasta.list \
  --output beta_results_30 \
  --cpu 35

#Annotate Family Representative Sequences via eggnog-mapper
ppanggolin fasta \
  -p beta_results_30/pangenome.h5 \
  --output beta_results_30/extracted_fasta \
  --prot_families all

#Map with COG annotations
mkdir beta_results_30/extracted_annotations 
emapper.py \
  -i beta_results_30/extracted_fasta/all_protein_families.faa \
  -o eggnog_annotations \
  --output_dir beta_results_30/extracted_annotations \
  --data_dir /RAID1/working/R324/Chiao/pan/eggnog_db \
  -m diamond \
  --pfam_realign realign \
  --cpu 35

##export gene family table
python eggnog_gene_family_table
#!/usr/bin/env python3
import re
from collections import defaultdict
import pandas as pd

# File paths
EGGNOG_FILE = "beta_results_30/extracted_annotations/eggnog_annotations.emapper.annotations"
GENE_FAMILIES_OUT = "beta_results_30/eggnog_gene_families.tsv"
FAMILY_MATRIX_OUT = "beta_results_30/eggnog_family_presence_absence.tsv"

print("Parsing eggNOG annotations into Gene Families...")

records = []
family_counts = defaultdict(lambda: defaultdict(int))

with open(EGGNOG_FILE, "r") as f:
    for line in f:
        # Skip header/comment lines
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        # Ensure standard eggNOG annotation row length
        if len(parts) < 12:
            continue

        query_id = parts[0].strip()
        ogs_str = parts[4].strip()         # Column 5: eggNOG_OGs
        cog_cats = parts[6].strip()        # Column 7: COG_category
        gene_name = parts[8].strip() if len(parts) > 8 and parts[8].strip() != "" else "-"
        description = parts[7].strip() if len(parts) > 7 and parts[7].strip() != "" else "-"
        pfams = parts[11].strip() if len(parts) > 11 and parts[11].strip() != "" else "-"

        # Skip unannotated rows
        if ogs_str == "-" or not ogs_str:
            continue

        # Extract Clean Genome ID (looks for GCA_xxxxxx.x or GCF_xxxxxx.x, or strips trailing CDS/gene IDs)
        gca_match = re.search(r"(GC[AF]_\d+\.\d+)", query_id)
        if gca_match:
            genome_id = gca_match.group(1)
        else:
            # Fallback regex if genome IDs follow non-NCBI naming conventions
            genome_id = re.sub(r"(_CDS_|_gene_|_pep_|_prot_|\.g\d+|\.\d+$).*", "", query_id)

        # Extract the finest Orthologous Group ID (last OG before @)
        og_list = ogs_str.split(",")
        finest_og = og_list[-1].split("@")[0]  # Extracts e.g. COG0457 or 1X5K2

        records.append({
            "gene_id": query_id,
            "genome_id": genome_id,
            "gene_family": finest_og,
            "full_OG_hierarchy": ogs_str,
            "COG_category": cog_cats if cog_cats != "" else "-",
            "gene_name": gene_name if gene_name != "" else "-",
            "Description": description if description != "" else "-",
            "PFAMs": pfams if pfams != "" else "-"
        })

        family_counts[finest_og][genome_id] += 1

# 1. Export Gene-to-Family Mapping Table with exact 8 headers
df_genes = pd.DataFrame(records)

# Ensure correct column ordering matching target table
df_genes = df_genes[[
    "gene_id", 
    "genome_id", 
    "gene_family", 
    "full_OG_hierarchy", 
    "COG_category", 
    "gene_name", 
    "Description", 
    "PFAMs"
]]

df_genes.to_csv(GENE_FAMILIES_OUT, sep="\t", index=False)

# 2. Export Gene Family Presence/Absence Matrix
df_matrix = pd.DataFrame.from_dict(family_counts, orient="index").fillna(0).astype(int)
df_matrix.to_csv(FAMILY_MATRIX_OUT, sep="\t")

print(f"\nSuccessfully generated Gene Family Tables!")
print(f" Total annotated genes mapped to families: {len(df_genes)}")
print(f" Total unique eggNOG gene families: {len(df_matrix)}")
print(f" Output 1 (Gene-level table): {GENE_FAMILIES_OUT}")
print(f" Output 2 (Family P/A matrix): {FAMILY_MATRIX_OUT}")

##map with subtilase and the syncom name
python summarize_subtilse_beta30.py

#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# File paths
INPUT_FILE = "beta_results_30/eggnog_gene_families.tsv"
OUTPUT_CSV = "beta_results_30/subtilase_summary_by_strain.csv"
OUTPUT_PLOT = "beta_results_30/subtilase_count_per_strain.png"

# Updated Strain Mapping Table (30 Strains)
GENOME_TO_STRAIN = {
    "GCA_001428165.1": "Root170",
    "GCA_001427425.1": "Root565",
    "GCA_001428845.1": "Root83",
    "GCA_001428645.1": "Root217",
    "GCA_001428665.1": "Root219",
    "GCA_001429215.1": "Root267",
    "GCA_001428965.1": "Root275",
    "GCA_001425125.1": "Root402",
    "GCA_001427445.1": "Root568",
    "GCA_001428225.1": "Root70",
    "GCA_001427305.1": "Root1480D1",
    "GCA_001428585.1": "Root198D2",
    "GCA_001425005.1": "Root336D2",
    "GCA_001428625.1": "Root209",
    "GCA_001426525.1": "Root133",
    "GCA_001427325.1": "Root1485",
    "GCA_001425685.1": "Root335",
    "GCA_001425045.1": "Root351",
    "GCA_001425265.1": "Root418",
    "GCA_001425905.1": "Root1272",
    "GCA_001428545.1": "Root189",
    "GCA_001425705.1": "Root1217",
    "GCA_001425055.1": "Root1237",
    "GCA_001426885.1": "Root1444",
    "GCA_001425165.1": "Root405",
    "GCA_001427705.1": "Root662",
    "GCA_001424835.1": "Root318D1",
    "GCA_001425205.1": "Root411",
    "GCA_001426595.1": "Root434",
    "GCA_001426505.1": "Root473",
}

print(f"Loading {INPUT_FILE}...")
df = pd.read_csv(INPUT_FILE, sep="\t")

# Map Accessions to Strain Names
df["strain_name"] = (
    df["genome_id"].map(GENOME_TO_STRAIN).fillna(df["genome_id"])
)

# Convert metadata columns to string for safe searching
df["full_OG_hierarchy"] = df["full_OG_hierarchy"].fillna("").astype(str)
df["gene_name"] = df["gene_name"].fillna("").astype(str)

# Search pattern targeting subtilase/subtilisin COG IDs and common gene symbols
SUBTILASE_PATTERN = r"COG1404|COG4934|subtilase|subtilisin|aprE|vpr|ispW|epiP|prtP"

# Filter rows matching subtilase terms
df_subtilase = df[
    df["full_OG_hierarchy"].str.contains(SUBTILASE_PATTERN, case=False, regex=True) |
    df["gene_name"].str.contains(SUBTILASE_PATTERN, case=False, regex=True)
].copy()

# 1. Total Subtilase Gene Count per Strain
summary = (
    df_subtilase.groupby("strain_name")
    .agg(
        total_subtilase_genes=("gene_id", "count"),
        unique_subtilase_families=("gene_family", "nunique"),
    )
    .reset_index()
    .sort_values(by="total_subtilase_genes", ascending=False)
)

print("\n--- Subtilase Family Breakdown per Strain ---")
print(summary.to_string(index=False))

# Save Summary Table to CSV
summary.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved summary table to: '{OUTPUT_CSV}'")

# 2. Plot Bar Chart of Subtilase Counts
if not summary.empty:
    # Scale width dynamically according to the number of strains present
    plot_width = max(12, len(summary) * 0.45)
    plt.figure(figsize=(plot_width, 6))
    
    ax = sns.barplot(
        data=summary,
        x="strain_name",
        y="total_subtilase_genes",
        palette="Oranges_r",
        edgecolor="black",
    )

    plt.title(
        "Total Subtilase (COG1404 / COG4934) Genes per Strain",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    plt.xlabel("Strain Name", fontsize=12, fontweight="bold")
    plt.ylabel("Gene Count", fontsize=12, fontweight="bold")
    
    # Rotate strain labels 45 degrees for better readability
    plt.xticks(rotation=45, ha="right", fontsize=10)

    # Add count labels on top of each bar
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(
                f"{int(height)}",
                (p.get_x() + p.get_width() / 2.0, height),
                ha="center",
                va="center",
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=9,
                fontweight="bold",
            )

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    plt.close()
    print(f"Saved bar plot to: '{OUTPUT_PLOT}'")
else:
    print("\nNo subtilase matches found. Check if COG1404 or COG4934 exist in your dataset.")






#for gamma
#Annotate Family Representative Sequences via eggnog-mapper
ppanggolin fasta \
  -p gamma_results_28/pangenome.h5 \
  --output gamma_results_28/extracted_fasta \
  --prot_families all

#Map with COG annotations
mkdir gamma_results_28/extracted_annotations 
emapper.py \
  -i gamma_results_28/extracted_fasta/all_protein_families.faa \
  -o eggnog_annotations \
  --output_dir gamma_results_28/extracted_annotations \
  --data_dir /RAID1/working/R324/Chiao/pan/eggnog_db \
  -m diamond \
  --pfam_realign realign \
  --cpu 35


