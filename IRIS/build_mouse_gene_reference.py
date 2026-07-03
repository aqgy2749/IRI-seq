from __future__ import annotations

import collections
import pickle
from pathlib import Path

import HTSeq
import pandas as pd


GTF = Path("/p2/zulab/jtian/data/IRISeq/reference/processed/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf")
OUT = Path("/p2/zulab/jtian/data/IRISeq/reference/Gene_annotation/mm10_GENCODE_M25_Gene_reference.pickle")


def attr(feature: HTSeq.GenomicFeature, key: str, default: str = "NA") -> str:
    return feature.attr.get(key, default)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    exon_only = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    gene_start: dict[str, set[int]] = collections.defaultdict(set)
    gene_end: dict[str, set[int]] = collections.defaultdict(set)

    gene_rows: dict[str, dict[str, object]] = {}
    exon_rows: list[dict[str, object]] = []
    gene_index = 0

    reader = HTSeq.GFF_Reader(str(GTF), end_included=True)
    for feature in reader:
        gene_id = attr(feature, "gene_id", "")
        if not gene_id:
            continue

        gene_type = attr(feature, "gene_type", attr(feature, "gene_biotype", "NA"))
        gene_name = attr(feature, "gene_name", gene_id)

        if feature.type == "gene":
            genes[feature.iv] += gene_id
            gene_start[gene_id].add(feature.iv.start_d)
            gene_end[gene_id].add(feature.iv.end_d)
            if gene_id not in gene_rows:
                gene_index += 1
                gene_rows[gene_id] = {
                    "Gene_id": gene_id,
                    "Gene_type": gene_type,
                    "Gene": "gene",
                    "Gene_name": gene_name,
                    "Index_ID": gene_index,
                }

        elif feature.type == "transcript":
            gene_start[gene_id].add(feature.iv.start_d)
            gene_end[gene_id].add(feature.iv.end_d)

        elif feature.type == "exon":
            exon_id = attr(feature, "exon_id", "No_exon_id")
            exons[feature.iv] += gene_id
            exon_only[feature.iv] += exon_id
            exon_rows.append(
                {
                    "Gene_id": gene_id,
                    "Gene_type": gene_type,
                    "Exon_id": exon_id,
                    "Gene_name": gene_name,
                }
            )

    gene_annotat = pd.DataFrame(gene_rows.values())
    exon_annotat = pd.DataFrame(exon_rows).drop_duplicates()

    gene_reference = {
        "genes": genes,
        "exons": exons,
        "exon_only": exon_only,
        "gene_end": dict(gene_end),
        "gene_start": dict(gene_start),
        "gene_annotat": gene_annotat,
        "exon_annotat": exon_annotat,
    }

    with OUT.open("wb") as handle:
        pickle.dump(gene_reference, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Wrote: {OUT}")
    print(f"Genes: {len(gene_annotat)}")
    print(f"Exons: {len(exon_annotat)}")


if __name__ == "__main__":
    main()
