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

  
