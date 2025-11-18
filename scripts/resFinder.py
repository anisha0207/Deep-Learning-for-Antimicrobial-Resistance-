import os
import subprocess
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from io import StringIO

FASTA_DIR = "fasta_files"
GBFF_DIR = "ncbi_dataset"
TSV_OUTPUT = "abricate_output.tsv"
CSV_OUTPUT = "resistance_summary.csv"

print("--> Running Abricate on .fna files...")

# Find all fasta files
fna_files = [f for f in os.listdir(FASTA_DIR) if f.endswith(".fna")]
fna_paths = [os.path.join(FASTA_DIR, f) for f in fna_files]

# Run abricate and save output
with open(TSV_OUTPUT, "w") as out_f:
    for fna in tqdm(fna_paths, desc="--> Scanning with Abricate"):
        result = subprocess.run(
            ["abricate", "--db", "resfinder", fna],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        out_f.write(result.stdout)

print("--> Parsing abricate TSV...")

# Extract header + data from TSV
with open(TSV_OUTPUT, "r") as f:
    lines = f.readlines()

header_line = [line for line in lines if line.startswith("#")][0].lstrip("#").strip()
data_lines = [line for line in lines if not line.startswith("#") and line.strip() != ""]

if not data_lines:
    raise ValueError("--> Abricate output missing required data rows!")

tsv_str = header_line + "\n" + "".join(data_lines)
abricate_df = pd.read_csv(StringIO(tsv_str), sep="\t")

required_cols = {"FILE", "SEQUENCE", "GENE", "PRODUCT"}
if not required_cols.issubset(abricate_df.columns):
    raise ValueError("--> Abricate output missing required columns!")

print("--> Matching GBFF gene info...")

# Group genes by file and SEQUENCE
res_lookup = abricate_df.groupby(["FILE", "SEQUENCE"])

results = []

for assembly_path, group in tqdm(res_lookup, desc="--> Processing genomes"):
    file_path, seq_id = assembly_path
    assembly_id = os.path.basename(file_path).replace(".fna", "")
    gbff_path = os.path.join(GBFF_DIR, assembly_id, "genomic.gbff")

    if not os.path.exists(gbff_path):
        print(f"--> Missing GBFF file for {assembly_id}")
        continue

    # Load full genome sequence
    with open(file_path, "r") as f:
        genome_seq = "".join(str(record.seq).upper() for record in SeqIO.parse(f, "fasta"))

    gc_count = genome_seq.count("G") + genome_seq.count("C")
    total_length = len(genome_seq)
    gc_percent = round(100 * gc_count / total_length, 2) if total_length else 0

    # Parse GBFF and find matching genes
    for record in SeqIO.parse(gbff_path, "genbank"):
        if record.id != seq_id:
            continue

        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers
            gene = qualifiers.get("gene", [""])[0]
            product = qualifiers.get("product", [""])[0]

            if not any(gene == row["GENE"] or gene in row["PRODUCT"] for _, row in group.iterrows()):
                continue

            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = feature.location.strand

            cds_seq = record.seq[start-1:end]
            if strand == -1:
                cds_seq = cds_seq.reverse_complement()

            start_codon = str(cds_seq[:3])
            end_codon = str(cds_seq[-3:])

            results.append({
                "Assembly_ID": assembly_id,
                "Seq_ID": record.id,
                "Gene": gene,
                "Product": product,
                "Start": start,
                "End": end,
                "GC": gc_count,
                "GC%": gc_percent,
                "Length": total_length,
                "Start_Codon": start_codon,
                "End_Codon": end_codon
            })

if not results:
    print("--> No resistance genes found.")
else:
    pd.DataFrame(results).drop_duplicates().to_csv(CSV_OUTPUT, index=False)
    print(f"--> Done! Results saved to: {CSV_OUTPUT}")
