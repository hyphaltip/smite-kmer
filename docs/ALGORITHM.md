# Algorithm: TE k-mer copy-number screen

This document describes what the code in `scripts/` does. It also states what has and has not been tested.

## Goal

Find transposable elements (TEs) that have more copies in one strain than in another. The input is short-read Illumina data for each strain and a FASTA file of TE sequences (`RiceTE.fa`). No genome assembly is needed.

## Idea

A k-mer from a TE appears in the reads once for every TE copy that carries it. A k-mer from a single-copy genome region appears once per genome. Both are sampled at the same read depth. So the depth of a TE k-mer, divided by the depth of a single-copy k-mer, estimates the TE copy number.

## Steps

### 1. Count k-mers (`scripts/count.sh`)

- One jellyfish index per strain. All R1 and R2 files of the strain go into the same index.
- k-mers are canonical (`-C`), so a k-mer and its reverse complement count together.
- Base quality: Phred+33, minimum Q20. A k-mer that contains a base below Q20 is not counted.
- If the sample sheet gives a folder, the script finds every read pair in it. R1 files match `*pair1*`, `*_R1*` or `*_1.f*`. The mate is found by name swap. A missing mate stops the run. Symlinks are followed.
- `jellyfish histo -h 1000` writes the k-mer count histogram.
- Run for k=21. The method also supports k=31. Only k=21 has been run.

### 2. Make the TE k-mer table (`te_kmer_table`)

- Take every k-mer of every TE sequence. Skip k-mers that contain a base other than A, C, G or T.
- Convert each to its canonical form.
- `multiplicity` is the number of times the canonical k-mer occurs in the whole TE file. This includes repeats inside one TE, such as terminal inverted repeats, and k-mers shared by two TEs.
- Only k-mers with `multiplicity == 1` are used for copy-number estimates. The full table, with multiplicities, is kept in the output.

### 3. Look up TE k-mer depths (`query_depths`)

`jellyfish query` gives the count of each canonical TE k-mer in the strain index. A k-mer that is absent has depth 0.

### 4. Find the single-copy depth (`baseline_peak`)

1. Smooth the histogram with a 3-bin moving average.
2. Start at count 2. Move up while the smoothed value keeps falling. This finds the end of the error valley.
3. The baseline is the count of the highest smoothed value at or above the valley.
4. If the peak is within 1 of the valley, the code prints a warning.

The rule takes the highest peak above the valley. It assumes that peak is the single-copy peak. It does not check this against a genome size or a heterozygous peak.

### 5. Estimate copy number (`analyze`)

For each strain and each unique TE k-mer:

    cn = depth / baseline

For each TE and strain, the copy number is the median of `cn` over the unique k-mers of that TE. Both depths are sampled at the same read depth, so the ratio does not need a correction for read length or k.

### 6. Compare strains to a reference

For each TE and each non-reference strain:

- Fold change `fc` = median `cn` of the strain / median `cn` of the reference.
- Paired Wilcoxon test over the unique k-mers (strain `cn` against reference `cn`).
- A 95% confidence interval for `fc` from a block bootstrap. The bootstrap picks blocks of k adjacent k-mer positions, because overlapping k-mers are not independent. It uses 2000 resamples and seed 1. It takes the 2.5th and 97.5th percentiles of the ratio of medians.
- The TE is flagged (`spike_<strain>`) when `fc` is at least 2.0 and the lower end of the confidence interval is above 1.

The Wilcoxon p-value is reported but is not used for the flag. There is no multiple-testing correction. Neighbouring k-mers overlap, so the Wilcoxon p-values are too small and should be read as a rough guide.

### 7. Position profile (`profile`)

This step shows where along a TE the depth is high or low.

- Split the unique k-mers of a TE into windows of 40 consecutive k-mers. A last window with fewer than 20 k-mers merges into the previous one.
- Window copy number = median depth in the window / baseline.
- `breadth` = fraction of windows with copy number at least 0.5 times the median window copy number.
- `min_window_ratio` = lowest window copy number / median window copy number.

`breadth` is meant to show partial-length copies. It is not calibrated. In the 12-strain run it was below 0.9 for 122 of 148 strain-TE pairs with median copy number of 2 or more. Do not use it as a flag yet.

### 8. GC-bias curve and correction (`gcbias`, `correct_gc`)

Sequencing depth can depend on the GC content of a k-mer. This step estimates that effect for each strain.

1. Dump the k-mers in the index whose count is between 0.5 and 1.5 times the baseline. These are mostly single-copy genome k-mers.
2. Group them by their number of G and C bases (0 to k).
3. For each group, `bias` = mean depth in the group / mean depth over all groups.
4. For a TE k-mer, corrected depth = depth / `bias` of its GC group. Groups with fewer than 1,000 genome k-mers get `bias` = 1.

`profile --gc-bias-dir` uses the corrected depths. The correction is not yet used in `analyze`. The correction has not yet been run on the real strains.

## Outputs

| File | Content |
|---|---|
| `baseline.k<k>.tsv` | single-copy depth per strain |
| `te_kmer_depths.k<k>.tsv` | every TE k-mer, with multiplicity and depth per strain |
| `te_copy_number.k<k>.tsv` | per TE: copy number per strain, `fc`, Wilcoxon p, bootstrap interval, spike flag |
| `te_tracks.k<k>.png` | copy number along each TE (15-k-mer rolling median) |
| `te_windows.k<k>.tsv`, `te_breadth.k<k>.tsv` | `profile` output |
| `<strain>.k<k>.gcbias.tsv` | GC-bias curve |

## What has been tested

- Unit tests cover canonical k-mers, multiplicity, the baseline peak, window profiles, GC counting and GC correction. A synthetic-read test with known copy numbers is in `tests/` and needs jellyfish.
- On 10 low-mPing rice accessions (Chen et al. 2019, Supplementary Data 1), the mPing copy number from this method correlates with the paper's read-depth estimate: Pearson r = 0.958, Spearman 0.957, n = 10. Our values are higher by a median of 1.17 times (range 0.94 to 1.48).
- Our NB estimate for mPing is 52.4. The paper lists 52 for its Nipponbare entry (CX140).
- Our HEG4 estimate is 485 copies. The paper reports 231 to 503 for the bursting accessions. It has no per-line value for HEG4 in the table used here. This number has no independent check.

## Known limits

- **Diverged copies.** A k-mer must match exactly. Copies with mutations lose k-mers. The method then underestimates copy number for old or diverged copies.
- **Shared k-mers.** `multiplicity` is counted only inside `RiceTE.fa`. A TE k-mer can also occur in other genome sequence, such as a related TE family not in the file. This inflates depth. It has not been measured here.
- **Baseline.** The baseline is the highest peak above the error valley. It is wrong if that peak is not the single-copy peak. The HEG4 baseline is 133 and NB is 57. HEG4 reads come from a folder named `HEG4_2-5kb`, which suggests a mate-pair library. This has not been checked, so the HEG4 baseline has not been confirmed.
- **Low depth.** For strains with a single-copy depth near 5 to 9, the error valley and the single-copy peak are close together. The baseline can be off by one or two counts, which is a large relative error.
- **Reference choice.** NB carries about 52 mPing copies. Fold change against NB understates how low other strains are. Absolute copy number, or a comparison to a set of low-copy strains, is more informative. A screen based on a low-copy set is planned and is not built.
- **GC dip in mPing.** Several strains show lower depth at mPing positions about 200 to 240, which have the highest GC (0.68). The cause is not yet known. The GC correction exists to test this.
