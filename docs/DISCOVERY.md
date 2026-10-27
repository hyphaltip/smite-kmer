# Blind discovery of TE bursts

Written 2026-10-27. This describes `scripts/discover_bursts.py` and what it has shown so far.

## Goal

Find sequences that are high copy in one strain and not in others, without a TE library. Rebuild the
sequence from the k-mers. mPing is the test case because its answer is known.

## Method

1. Dump the test strain's k-mers at or above `--min-cn` (20) times its single-copy depth.
2. Look each up in the comparison strains. Keep k-mers whose copy number is at least `--min-fold` times the
   comparison median plus 1. Optional hysteresis (`--extend-fold`) keeps neighbours at a lower fold if the
   unitig has enough strict seed k-mers.
3. Build a de Bruijn graph (odd k, canonical k-mers). Drop edges to k-mers rarer than 0.1 of their siblings
   (sequencing errors). Pop bubbles (two paths that rejoin within 2k steps, as from a common second allele).
   Read off unitigs.
4. Report each unitig with its copy number in the test strain and in every comparison strain, its fold,
   `n_carriers` (comparison strains with at least 3 copies and a quarter of the test copy number), and a flag
   for sequences made mostly of masked k-mers (adapters, PhiX). Output is ranked by fold times length.

`RiceTE.fa` is used only to label results. It does not affect which sequences are found.
The test strain needs a single-copy depth of at least 15 (`--min-depth`). Below that, baselines are unreliable.

## Results (k=21)

Test strains: HEG4, A119, A123, A157 (burst lines), and B001, B269, IRIS_313-11689, NB (52 to 71 mPing copies in
Chen et al. 2019). Negative control: IRIS_313-11652 (5 copies). These are the only strains with an index and a
single-copy depth of at least 15.

- mPing was recovered as one unitig in every test strain. Where the sequence was compared with the library,
  all 430 bp match with 0 mismatches. The unitigs are 432 to 434 bp.
- The negative control gave no mPing.
- With 10 comparison strains, the non-mPing, non-artifact unitigs per test strain were 0 to 32. With 34
  comparison strains across six subgroups, they were 0 to 24 in the runs without hysteresis. More of them then had
  carriers among the comparison strains (IRIS_313-11652: 12 of 18 private unitigs became 1 of 13).
- Hysteresis did not help on these data: unitig counts rose (HEG4 24 to 32, NB 24 to 40) and mPing was
  already complete without it.
- No unitig matched the rice plastid (NC_001320.1) or mitochondrion (NC_011033.1).

### What the other unitigs are (carrier survey over 80 indexed strains, `scripts/carrier_survey.py`)

- Illumina adapter (flagged). Present in the Wessler-lab libraries (HEG4, NB).
- Two ~340 bp unitigs (GC 0.68 to 0.70) in HEG4 and NB, present at 5 or more copies in 64 to 69 of 80 strains of
  all subgroups, at 240 to 280 copies in HEG4. A reviewer found 5.8S rRNA-like sequence in them. I have not
  confirmed this. They look like rDNA, a high-copy family with variable copy number, not a recent burst.
- Three unitigs (88 to 136 bp) in HEG4 only, at 150 to 419 copies, and at about 10 copies in A119 and A123. They
  contain cloning polylinker sites (GGATCC CCCGGG, GCGGCCGC GATATC CTGCAG). They match neither the Illumina nor the
  Nextera mate-pair adapter files. This looks like vector or linker DNA. It has not been checked against a vector
  database.
- IRIS_313-11652: six GC-rich unitigs at 800 to 1,280 copies. Six strains carry them (five temperate japonica,
  one admixed). This looks like a repeat shared within one lineage, not contamination. This is an inference.

## Limits

- Not tested: whether a burst of a different kind (long LTR element, diverged family) is found. mPing has
  near-identical copies, which suits exact k-mers.
- Lineage-specific repeats pass as "bursts". The comparison set must span the subgroups.
- The fold cutoff is fixed. mPing at 52 to 71 copies fails a cutoff of 10 because the comparison median is about
  6.5 copies.
- Graphs are fragmented. One run had 609 single-k-mer unitigs out of 985, because there is no tip removal.
- Only 9 strains can be test strains. The others are below 15x single-copy depth.
- The technical-replicate control (CX140 against NB) could not run, because CX140's index was not kept.
- Speed: HEG4 against 10 comparison strains took 260 s. Against 34 strains, with 8 parallel lookups, it took 446 s.
  Most of the time is `jellyfish query`. The graph step takes seconds.

## Not built

- Read recruitment, local assembly, flanks, target-site duplications and insertion sites.
- An empirical (MAD-based) threshold in place of the fixed fold.
- A cohort-wide k-mer by strain matrix, and a smaller persistent table per strain.
- A spike-in test of detection limits.

## Files

- `scripts/discover_bursts.py`, `scripts/discover_task.sh`, `scripts/carrier_survey.py`
- `tests/test_discover_bursts.py`
- `resources/masks/` (plastid and mitochondrion FASTA from NCBI)
- Results (not in git): `results/discovery_v2/`, `results/discovery_v3/`, `work/carrier_survey.tsv`
