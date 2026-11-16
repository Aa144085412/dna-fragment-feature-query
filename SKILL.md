---
name: dna-fragment-feature-query
description: Locate mutation-labelled SnapGene DNA fragments on human hg38 with local BLAST+, retrieve UCSC RepeatMasker/cCRE/RefSeq features, and produce evidence-ranked repair recommendations from SnapGene and optional vendor XLSX. Use for synthetic DNA fragment discrepancy reviews, mutation mapping and genomic feature queries.
compatibility: Python 3.10+, local NCBI BLAST+ and indexed hg38; UCSC API access only with opt-in
metadata:
  version: "0.1.0"
---

# DNA Fragment Feature Query

Invoke the **existing** BLAST+ and hg38 database (such as those already managed by SeqAtlas/PrimerForge); do not deploy a second indexing service or new Agent controller.

## Supported task

Read mutation-labelled `misc_feature` entries from one or more SnapGene `.dna` files, optionally reconcile vendor Excel (column B fragment IDs, E mutation descriptions), map ±150 bp to **hg38**, query UCSC `rmsk`, `cCREregistry`, `ncbiRefSeq`, assess reported coding SNP changes and annotate each mutation with one of the bounded repair recommendations.

## Inputs and execution

- Directory of `.dna` files, optionally matching annotation `.xlsx`.
- Known local hg38 BLAST database prefix and NCBI `blastn` binary.
- Explicit authorization for `--online` UCSC coordinate queries. Without it, report `LOCAL_MAPPING_ONLY` and **do not infer functional absence**.

Run, from this directory:

```bash
python scripts/dfq.py --input INPUT_FOLDER --xlsx VENDOR.xlsx \
  --blastn BLASTN_PATH --db HG38_BLAST_PREFIX --online --out OUTPUT_FOLDER
```

`--xlsx` is optional. See README for exact Windows paths.

## Mandatory scientific gates

1. **Separate report from observation.** A feature saying `1-A-G` is a claim from the vendor/annotation, not a confirmed sequencing genotype. `.dna` may contain the reference base. Do not call observed mutations or base editing without raw sequencing evidence.
2. Require sufficiently long, high-identity BLAST alignment and assess alternate contigs and ties. A multimap is not a unique mapping.
3. Output coordinates **1-based**; UCSC API intervals are **0-based half-open**.
4. Query hg38 reference nucleotide and UCSC **rmsk / ENCODE4 cCRE registry / RefSeq** with explicit success/failure for each site.
5. A gene-body overlap does not prove an exon/coding change. Apply CDS/exon/strand and transcript-specific codon calculations for predicted coding substitutions.
6. Repeat indels require trace review and normalization; cCRE is candidate evidence not proven enhancer function. Missing tracks, mapping uncertainty or reference-allele conflict must remain unresolved.
7. Never alter user DNA files; never push user sample data, confidential Excel, genome index, credentials or local result files to GitHub.

## Deliverables

`report.md`, `review.xlsx`, `review.tsv`, `result.json` with per-site BLAST alternatives, coordinates, annotations, confidence/uncertainty, exact reason and UCSC links. If `--xlsx` was supplied, **after** these reports, generate `<original-name>_hg38修复说明.xlsx`: preserve the original workbook structure and other sheets, fill F (repair / do-not-repair / verify) and G (per-site reasons), and never overwrite the input or an existing output. Uncertain evidence must remain `待核实` rather than become `不修复`.

Accept completion only if input counts reconcile, all site coordinates are traceable, API status is explicit, and outputs have been validated. If inputs/data are insufficient, state the smallest missing evidence and stop rather than invent a conclusion.
