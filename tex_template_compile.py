from datetime import date
from pathlib import Path
import argparse
import csv
import json

# Function to process CNV file
def process_cnv_file(file_path):
    cnv_data = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            (
                chr_info,
                numsnp_info,
                length_info,
                state_info,
                file_info,
                startsnp_info,
                endsnp_info,
                conf,
                iscn,
            ) = parts
            region = chr_info.split(":")
            chr_start_end = (
                f"{region[0]}_{region[1].split('-')[0]}_{region[1].split('-')[1]}"
            )
            cnv = {
                "iscn": iscn,
                "chr_info": chr_info,
                "numsnp_info": numsnp_info,
                "length_info": length_info,
                "state_info": state_info,
                "file_info": file_info,
                "startsnp_info": startsnp_info,
                "endsnp_info": endsnp_info,
                "conf": conf,
            }
            cnv_data[chr_start_end] = cnv
    return cnv_data

# Function to process Scoresheet file
def process_scoresheet_file(file_path):
    scoresheet_data = {}
    with open(file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variant_id = row["VariantID"]
            # Remove _DEL and _DUP from VariantID
            clean_variant_id = variant_id.replace("_DEL", "").replace("_DUP", "")
            scoresheet_data[clean_variant_id] = dict(row)
    return scoresheet_data

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Join CNV and Scoresheet files.")
    parser.add_argument("--cnv_file", help="Path to the CNV file")
    parser.add_argument("--scoresheet_file", help="Path to the Scoresheet file")
    parser.add_argument("--tex_template", help="Path to the tex file")
    parser.add_argument("--plot_dir", help="Path to the plot directory")
    parser.add_argument("--protocol_id", help="protocol id")
    parser.add_argument("--institute", help="institute")
    parser.add_argument("--chip_type", help="chip type", default="GSA-Cyto")
    parser.add_argument("--version", help="version", default="24.12")
    parser.add_argument("--sample_id", help="sample_id")
    parser.add_argument("--output_file", help="Path to the output tex file")
    parser.add_argument("--cnvs", help="Path to the premerged CNV file with scoresheet")

    args = parser.parse_args()

    if args.cnvs:
        chip_id, position = args.sample_id.split("_")
        with open(args.cnvs) as f:
            cnvs = json.load(f)
    else:
        cnv_path = Path(args.cnv_file)
        chip_id = cnv_path.stem.split("_")[0]
        position = cnv_path.stem.split("_")[1]
        # Read and process the CNV file
        cnv_data = process_cnv_file(args.cnv_file)

        # Read and process the Scoresheet file
        scoresheet_data = process_scoresheet_file(args.scoresheet_file)

        cnvs = {}
        for variant_id, cnv_dict in cnv_data.items():
            score_dict = scoresheet_data.get(variant_id, {})
            merged_dict = {**cnv_dict, **score_dict}
            cnvs[variant_id] = merged_dict

    latex_string = ""
    # Özet tablosunu yap
    for variant_id, cnv in cnvs.items():
        pass  # Placeholder for any summary table logic

    # Her varyant için ayrı ayrı subsectionlar
    for variant_id, cnv in cnvs.items():
        if "Chromosome" not in cnv.keys():
            continue

        num_snp = cnv["numsnp_info"].split("=")[1]
        cnv_length = cnv["length_info"].split("=")[1]
        cnv_conf = cnv["conf"].split("=")[1]
        iscn = cnv["iscn"].replace("_", "\\_")

        evidences = [
            "1A-B",
            "2A",
            "2B",
            "2C",
            "2D",
            "2E",
            "2F",
            "2G",
            "2H",
            "2I",
            "2J",
            "2K",
            "2L",
            "3",
            "4A",
            "4B",
            "4C",
            "4D",
            "4E",
            "4F-H",
            "4I",
            "4J",
            "4K",
            "4L",
            "4M",
            "4N",
            "4O",
            "5A",
            "5B",
            "5C",
            "5D",
            "5E",
            "5F",
            "5G",
            "5H",
        ]
        evidence_in_report = {}
        for evidence in evidences:
            evidence_value = cnv.get(evidence, 0)
            if evidence_value == '':
                evidence_value = 0.0
            if float(evidence_value) != 0:
                evidence_in_report[evidence] = evidence_value
        evidence_str = " ".join([f"{k}: {v}" for k, v in evidence_in_report.items()])

        # Handle gene names: italicize and handle empty cases
        known_genes = cnv.get('Known or predicted dosage-sensitive genes', '').strip()
        all_genes = cnv.get('All protein coding genes', '').strip()

        # Function to italicize gene names
        def italicize_genes(gene_str):
            if not gene_str:
                return "Gen içermemektedir."
            genes = [gene.strip() for gene in gene_str.split(",")]
            italic_genes = ", ".join([f"\\textit{{{gene}}}" for gene in genes if gene])
            return italic_genes if italic_genes else "Gen içermemektedir."

        known_genes_str = italicize_genes(known_genes)
        all_genes_str = italicize_genes(all_genes)

        latex_string += f"""
\\subsection{{{variant_id.replace('_', ':', 1).replace('_', '-')}}}

\\begin{{tabularx}}{{\\textwidth}}{{l X X X X X}}
ISCN           &  State         & SNP number   & Length        & ACMG classification     & Confidence    \\\\
\\hline
{iscn}         & {cnv["Type"]}  & {num_snp}    & {cnv_length}  & {cnv["Classification"]} & {cnv_conf}     \\\\
\\hline
\\end{{tabularx}}

\\vspace{{0.5cm}}  % Bu kısmı boşluk eklemek için kullanıyoruz

\\begin{{tabularx}}{{\\textwidth}}{{l L}}  % Tablo genişliği otomatik ayarlanır
ACMG Total Score:                                        & {cnv.get("Total score", "N/A")} \\\\
Evidence:                                                & {evidence_str} \\\\
Known or Predicted Dosage-Sensitive Genes:  & {known_genes_str} \\\\
All Protein Coding Genes:                   & {all_genes_str}
\\end{{tabularx}}

\\vspace{{0.5cm}}  % Bu kısmı boşluk eklemek için kullanıyoruz

\\begin{{center}}
\\includegraphics[width=0.8\\textwidth]{{ {args.plot_dir}/plot_{variant_id.replace("chr", "").replace("_DEL", "").replace("_DUP", "")}.png }}
\\end{{center}}
"""

    with open(args.tex_template, "r", encoding="utf-8") as in_file, open(
        args.output_file, "w", encoding="utf-8"
    ) as out_file:

        output_template = in_file.read().replace("%%BULGULAR%%", latex_string)
        output_template = output_template.replace("%%genome_plot%%", f"{args.plot_dir}/plot_genome.png")
        output_template = output_template.replace("%%institute%%", args.institute)
        output_template = output_template.replace("%%canvasVersion%%", args.version)
        output_template = output_template.replace("%%protocolId%%", args.protocol_id)
        output_template = output_template.replace(
            "%%summaryDate%%", date.today().strftime("%B %d, %Y")
        )
        output_template = output_template.replace("%%chipId%%", chip_id)
        output_template = output_template.replace("%%chipType%%", args.chip_type)
        output_template = output_template.replace("%%chipPosition%%", position)
        out_file.write(output_template)

if __name__ == "__main__":
    main()

