import os
from Bio import SeqIO
from tqdm import tqdm  # <-- progress bar

input_dir = "/Users/sanjoydasgupta/Desktop/genomics_env/dna_playground/scripts/ncbi_dataset"
output_dir = "fasta_files"
os.makedirs(output_dir, exist_ok=True)

# Collect all .gbff file paths
gbff_files = []
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if file.endswith(".gbff"):
            gbff_files.append(os.path.join(root, file))

print(f"🧬 Found {len(gbff_files)} .gbff files to convert\n")

# Process each .gbff file with progress bar
for gbff_path in tqdm(gbff_files, desc="Converting to FASTA", unit="file"):
    asm_id = os.path.basename(os.path.dirname(gbff_path))  # use parent folder name
    fasta_path = os.path.join(output_dir, f"{asm_id}.fna")
    
    with open(fasta_path, "w") as out_f:
        for record in SeqIO.parse(gbff_path, "genbank"):
            SeqIO.write(record, out_f, "fasta")

print("\n✅ All .gbff files successfully converted to .fna FASTA format.")

