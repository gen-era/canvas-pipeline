from datetime import date
from pathlib import Path
import argparse
import csv
import json
import re


def escape_latex(text):
    if not isinstance(text, str):
        return text

    latex_special_chars = {
        "&": "\\&",
        "%": "\\%",
        "$": "\\$",
        "#": "\\#",
        "_": "\\_",
        "{": "\\{",
        "}": "\\}",
        "~": "\\textasciitilde{}",
        "^": "\\textasciicircum{}",
        "\\": "\\textbackslash{}",
        "<": "\\textless{}",
        ">": "\\textgreater{}",
        "|": "\\textbar{}",
        '"': "\\textquotedbl{}",
        "`": "\\textasciigrave{}",
        "'": "\\textquotesingle{}",
    }

    text = text.replace("\\", "\\textbackslash{}")

    for char, escape in latex_special_chars.items():
        if char != "\\":
            text = text.replace(char, escape)

    return text


def get_confidence(cnv: dict) -> str:

    raw = None

    for k in ("conf", "Confidence", "confidence"):
        if k in cnv and cnv[k] not in (None, ""):
            raw = str(cnv[k])
            break

    if raw is None:
        return "-"

    if "=" in raw:
        raw = raw.split("=", 1)[1]

    val = raw.strip()

    if val.lower() in ("undefined", "na", "n/a", "none", ""):
        return "-"

    return val


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


def process_scoresheet_file(file_path):

    scoresheet_data = {}

    with open(file_path, "r") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            variant_id = row["VariantID"]

            clean_variant_id = variant_id.replace("_DEL", "").replace("_DUP", "")

            scoresheet_data[clean_variant_id] = dict(row)

    return scoresheet_data


def main():

    parser = argparse.ArgumentParser(description="Join CNV and Scoresheet files.")

    parser.add_argument("--cnv_file")
    parser.add_argument("--scoresheet_file")
    parser.add_argument("--tex_template")
    parser.add_argument("--plot_dir")
    parser.add_argument("--protocol_id")
    parser.add_argument("--institute")
    parser.add_argument("--chip_type", default="GSA-Cyto")
    parser.add_argument("--version", default="24.12")
    parser.add_argument("--sample_id")
    parser.add_argument("--output_file")
    parser.add_argument("--cnvs")

    parser.add_argument(
        "--include_plots",
        action="store_true",
        help="Include LogR/BAF plots in report",
    )

    args = parser.parse_args()

    chipsample_notes = None

    if args.cnvs:

        chip_id, position = args.sample_id.split("_")

        with open(args.cnvs) as f:

            cnvs = json.load(f)

            chipsample_notes = cnvs["chipsample_notes"]

            cnvs = cnvs["cnvs"]

    else:

        cnv_path = Path(args.cnv_file)

        chip_id = cnv_path.stem.split("_")[0]

        position = cnv_path.stem.split("_")[1]

        cnv_data = process_cnv_file(args.cnv_file)

        scoresheet_data = process_scoresheet_file(args.scoresheet_file)

        cnvs = {}

        for variant_id, cnv_dict in cnv_data.items():

            score_dict = scoresheet_data.get(variant_id, {})

            merged_dict = {**cnv_dict, **score_dict}

            cnvs[variant_id] = merged_dict

    latex_string = ""

    quality_notes = []
    bulgular_rows = []
    genes_blocks = []
    plot_blocks = []

    def italicize_genes(gene_str):

        if not gene_str:
            return "Gen içermemektedir."

        genes = [gene.strip() for gene in gene_str.split(",")]

        italic_genes = ", ".join(
            [f"\\textit{{{gene}}}" for gene in genes if gene]
        )

        return italic_genes if italic_genes else "Bölge Gen içermemektedir."

    if not cnvs:

        latex_string += "Yapılan mikroarray analizi sonucunda gözlenen kopya sayısı değişimleri arasında klinik önemi olan bir bulgu gözlenmemiştir."

    else:

        for variant_id, cnv in cnvs.items():

            if "Chromosome" not in cnv.keys():
                continue

            cnv_length = cnv["length_info"].split("=")[1]

            cnv_conf = get_confidence(cnv)

            iscn = cnv["iscn"].replace("_", "\\_")

            evidences = [
                "1A-B","2A","2B","2C","2D","2E","2F","2G","2H","2I","2J",
                "2K","2L","3","4A","4B","4C","4D","4E","4F-H","4I","4J",
                "4K","4L","4M","4N","4O","5A","5B","5C","5D","5E","5F","5G","5H"
            ]

            evidence_in_report = {}

            for evidence in evidences:

                evidence_value = cnv.get(evidence, 0)

                if evidence_value == "":
                    evidence_value = 0.0

                if float(evidence_value) != 0:
                    evidence_in_report[evidence] = evidence_value

            evidence_str = " ".join(
                [f"{k}: {v}" for k, v in evidence_in_report.items()]
            )

            known_genes = cnv.get(
                "Known or predicted dosage-sensitive genes", ""
            ).strip()

            all_genes = cnv.get(
                "All protein coding genes", ""
            ).strip()

            known_genes_str = italicize_genes(known_genes)

            all_genes_str = italicize_genes(all_genes)

            bulgu_note = ""

            if "notes" in cnv and cnv["notes"]:
                bulgu_note = " ".join(
                    escape_latex(note) for note in cnv["notes"]
                )

            type_safe = str(cnv.get("Type", "")).replace("_", "\\_")

            class_safe = str(cnv.get("Classification", "")).replace("_", "\\_")

            cnv_length_safe = str(cnv_length).replace("_", "\\_")

            bulgular_rows.append(
                (iscn, type_safe, cnv_length_safe, class_safe)
            )

            genes_blocks.append(
                (bulgu_note, known_genes_str, all_genes_str)
            )

            total_score = str(cnv.get("Total score", "N/A")).replace("_", "\\_")

            evidence_safe = (evidence_str or "—").replace("_", "\\_")

            cnv_conf_safe = str(cnv_conf).replace("_", "\\_")

            quality_notes.append(
                {
                    "iscn": iscn,
                    "evidences": evidence_safe,
                    "classification_score": f"{class_safe} ({total_score})",
                    "confidence": cnv_conf_safe,
                }
            )

            plot_blocks.append(
                f"""
\\vspace{{0.5cm}}
\\begin{{center}}
\\includegraphics[width=0.8\\textwidth]{{ {args.plot_dir}/plot_{variant_id.replace("chr","").replace("_DEL","").replace("_DUP","")}.png }}
\\end{{center}}
"""
            )

        latex_string += """
\\noindent
\\begin{tabularx}{\\textwidth}{m{1.4cm} X X >{\\centering\\arraybackslash}m{2.2cm} X}
\\textbf{Bulgu No} & \\textbf{ISCN} & \\textbf{Değişim Türü} & \\makecell{\\textbf{Değişim Boyutu}\\\\\\textbf{(bp)}} & \\textbf{Sınıflandırma} \\\\
\\hline
"""

        for no, (iscn, type_safe, cnv_length_safe, class_safe) in enumerate(
            bulgular_rows, start=1
        ):

            latex_string += (
                f"{no} & {iscn} & {type_safe} & {cnv_length_safe} & {class_safe} \\\\\n\\hline\n"
            )

        latex_string += "\\end{tabularx}\n\n"

        for no, (bulgu_note, known_genes_str, all_genes_str) in enumerate(
            genes_blocks, start=1
        ):

            latex_string += f"""
\\subsubsection{{Bulgu No {no}:}}
"""

            if bulgu_note:

                latex_string += f"""
\\noindent\\textbf{{Not:}} {bulgu_note}

"""

            latex_string += f"""
\\noindent\\textbf{{Bilinen veya tahmin edilen Dosage-Sensitive Genler:}}

{known_genes_str}

\\noindent\\textbf{{Protein kodlayan genler:}}

{all_genes_str}

\\vspace{{0.5cm}}
"""

        if args.include_plots:

            for plot_block in plot_blocks:

                latex_string += plot_block

        if quality_notes:

            latex_string += """
\\newpage
\\section{Kalite Notu}
"""

            latex_string += """
\\begin{tabularx}{\\textwidth}{m{1.4cm} X X >{\\centering\\arraybackslash}m{2.2cm}}
\\textbf{Bulgu No} & \\textbf{Kanıtlar} & \\textbf{Sınıflandırma (Skor)} & \\makecell{\\textbf{Kanıt}\\\\\\textbf{Güven Düzeyi}} \\\\
\\hline
"""

            for no, q in enumerate(quality_notes, start=1):

                latex_string += (
                    f"{no} & {q['evidences']} & {q['classification_score']} & {q['confidence']} \\\\\n\\hline\n"
                )

            latex_string += "\\end{tabularx}\n\n"

    with open(args.tex_template, "r", encoding="utf-8") as in_file, open(
        args.output_file, "w", encoding="utf-8"
    ) as out_file:

        latex_string = re.sub(r"(?<!\\)_", r"\\_", latex_string)

        output_template = in_file.read().replace(
            "%%BULGULAR%%", latex_string
        )

        if chipsample_notes:

            escaped_notes = "".join(
                escape_latex(note) for note in chipsample_notes
            )

            output_template = output_template.replace(
                "%%CHIPSAMPLENOTES%%",
                "\\textbf{{Not:}} " + escaped_notes,
            )

        else:

            output_template = output_template.replace(
                "%%CHIPSAMPLENOTES%%", ""
            )

        output_template = output_template.replace(
            "%%genome_plot%%",
            f"{args.plot_dir}/plot_genome.png",
        )

        output_template = output_template.replace(
            "%%institute%%",
            args.institute,
        )

        output_template = output_template.replace(
            "%%canvasVersion%%",
            args.version,
        )

        output_template = output_template.replace(
            "%%protocolId%%",
            args.protocol_id,
        )

        output_template = output_template.replace(
            "%%summaryDate%%",
            date.today().strftime("%B %d, %Y"),
        )

        output_template = output_template.replace(
            "%%chipId%%",
            chip_id,
        )

        output_template = output_template.replace(
            "%%chipType%%",
            args.chip_type,
        )

        output_template = output_template.replace(
            "%%chipPosition%%",
            position,
        )

        out_file.write(output_template)


if __name__ == "__main__":
    main()
