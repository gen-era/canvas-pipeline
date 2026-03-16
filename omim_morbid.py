import pandas as pd
from pathlib import Path
import functools

BASE = Path(__file__).resolve().parent


@functools.lru_cache(maxsize=1)
def load_tables():

    omim = pd.read_csv(BASE / "ucsc_omim_genes.tsv", sep="\t")

    gene_info = pd.read_csv(
        BASE / "Homo_sapiens.gene_info",
        sep="\t",
        low_memory=False
    )

    mim = pd.read_csv(
        BASE / "mim2gene_medgen",
        sep="\t",
        comment="#",
        header=None,
        names=[
            "MIM_number","GeneID","type",
            "Source","MedGenCUI","Comment"
        ]
    )

    def extract_mim(dbxrefs):
        if pd.isna(dbxrefs):
            return None
        for x in str(dbxrefs).split("|"):
            if x.startswith("MIM:"):
                return x.replace("MIM:", "")
        return None

    gene_info["omim_id"] = gene_info["dbXrefs"].apply(extract_mim)

    morbid = mim[
        (mim["type"] == "phenotype") &
        (mim["Source"].str.contains("GeneMap", na=False))
    ]

    return omim, gene_info, morbid


def get_omim_morbid_for_region(chr, start, end):

    omim, gene_info, morbid = load_tables()

    region = omim[
        (omim["chr"] == chr) &
        (omim["start"] <= end) &
        (omim["end"] >= start)
    ]

    if region.empty:
        return [], []

    ids = region["omim_id"].astype(str).unique()

    genes = gene_info[
        gene_info["omim_id"].isin(ids)
    ][["Symbol","GeneID"]]

    omim_genes = sorted(genes["Symbol"].dropna().unique())

    morbid_ids = morbid["GeneID"].astype(str)

    morbid_genes = sorted(
        genes[
            genes["GeneID"].astype(str).isin(morbid_ids)
        ]["Symbol"].unique()
    )

    return omim_genes, morbid_genes


if __name__ == "__main__":

    omim, morbid = get_omim_morbid_for_region(
        "chr5",
        31679,
        25936763
    )

    print("OMIM:", ", ".join(omim))
    print("Morbid:", ", ".join(morbid))
