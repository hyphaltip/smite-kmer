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

## Cohort-wide matrix (`scripts/cohort_matrix.py`)

One row per candidate k-mer and one column per strain, so a strain is judged against the whole cohort and
the carriers of a sequence can be counted. Candidates are the k-mers at or above 20 times the single-copy depth in
at least one strain: 2,339,175 k-mers across 80 strains (k=21). A k-mer that never reaches 20 times in any strain is
not a row, so families that are partly below that level are only partly present (for example 365 of 417 unique
Dasheng k-mers).

Stages: `strains`, `dump`, `union`, `query`, `build`, `call`, `evaluate`. Run through `scripts/cohort_task.sh`
(array stages in 8 parts) and `scripts/cohort_calls_task.sh`. The matrix is `work/matrix/matrix.npy` (float32, about
0.75 GB). Its mPing values equal the earlier analysis exactly (HEG4 484.6, NB 52.4, B031 12.4, A119 360.2).

Time on the stajichlab queue: dump 3 min, union 1 min, query 2.5 min, build 0.5 min (each of the array stages ran
as 8 parts). A `call` for one strain takes about 20 s once the matrix exists. The earlier design needed 4 to 14 min
per test strain.

### Results

- Detection limit (`evaluate`, rule: copy number >= 20 and >= 5 times the cohort median plus 1; the cohort
  median for mPing is 17.3): mPing stands out in 7 strains. They are the four lab lines (HEG4 485, A119 360,
  A123 276, A157 143 copies) and the three highest-copy 3K accessions (B160 paper value 180, B235 113, B005 86;
  91 to 98% of their mPing k-mers called). The strains with 54 to 71 paper copies are not called (8% of k-mers or
  fewer). The limit is about 90 copies. It follows from the cohort median, so a cohort with a lower median has a lower
  limit.
- This cohort is enriched for mid and high mPing because the 3K accessions were sampled across mPing bins. The
  paper table gives a median of 6 to 10 copies across the 3,000 lines. A cohort drawn at random would have a lower
  median and detect lower copy numbers.
- The carrier count does not separate a burst from a repeat shared by a lineage. mPing in A119 has 6 other
  carriers (the related lab lines and top 3K accessions), and the temperate japonica repeat in IRIS_313-11652 has 5. A
  limit of 5 carriers removed mPing from A119, A123 and A157 and left the lineage repeat. Burst lines that descend from
  one ancestor share the burst, so the carrier count is information to read, not a filter. Structure (terminal repeats,
  target-site duplications, a transposase) is still needed.
- At fold 10, non-mPing, non-artifact unitigs: HEG4 26 (10 comparison strains), 24 (34), 17 (matrix, 80);
  IRIS_313-11652 15 (10), 11 (matrix). mPing is complete in HEG4, A119 and A123 and missed at fold 10 in A157 (148
  copies against a cutoff of about 183).

## Not built

- Read recruitment, local assembly, flanks, target-site duplications and insertion sites.
- An empirical (MAD-based) threshold in place of the fixed fold.
- A smaller persistent per-strain table, so jellyfish indexes can be deleted.
- A spike-in test of detection limits.

## Files

- `scripts/discover_bursts.py`, `scripts/discover_task.sh`, `scripts/carrier_survey.py`
- `scripts/cohort_matrix.py`, `scripts/cohort_task.sh`, `scripts/cohort_calls_task.sh`
- `tests/test_discover_bursts.py`, `tests/test_cohort_matrix.py`
- `resources/masks/` (plastid and mitochondrion FASTA from NCBI)
- Results (not in git): `results/discovery_v2/`, `results/discovery_v3/`, `work/carrier_survey.tsv`
