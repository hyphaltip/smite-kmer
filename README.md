# smite-kmer

Screen short-read libraries for transposable-element (TE) copy-number differences
between strains using jellyfish k-mer counts. Built to ask whether mPing and other
elements in `RiceTE.fa` are amplified in one rice strain relative to another.

## Method

1. One jellyfish index per strain (R1 + R2 together), canonical k-mers, Phred+33
   quality mask (Q20). Run for k=21 and k=31.
2. Single-copy depth per strain = modal peak of its k-mer histogram above the error valley.
3. Canonical k-mers from `RiceTE.fa`; only k-mers occurring exactly once in the file are used
   for the primary estimate. The per-k-mer table keeps a multiplicity column.
4. Per-TE copy number = median(TE k-mer depth) / single-copy depth, per strain.
5. Fold-change vs a reference strain (`--reference`, default NB), a paired Wilcoxon test, and a
   block-bootstrap CI (blocks of k adjacent k-mers, since overlapping k-mers are not independent).
   A TE is flagged when fold-change >= 2 and the CI excludes 1.

## Usage

    mamba create -n mitekmer -c conda-forge python=3.12 numpy scipy pandas matplotlib pytest
    # jellyfish >= 2.3 on PATH (bioconda on Linux; `brew install jellyfish` on Apple silicon)

    # samples.tsv: strain <TAB> reads_dir   (all read pairs in the folder are combined)
    # or:           strain <TAB> r1 <TAB> r2 (a single pair)
    scripts/run_all.sh 21 8 2G NB      # K THREADS HASH_SIZE REFERENCE
    # on SLURM: one count job per strain, then an analysis job (afterok dependency)
    scripts/submit.sh 21 NB            # K REFERENCE; env: PARTITION ACCOUNT THREADS MEM TIME HASH_SIZE
    python -m pytest tests

In a reads folder, R1 files match `*pair1*`, `*_R1*` or `*_1.f*`. The mate is found by swapping
to `pair2`, `_R2` or `_2.`. A missing mate stops the run. All pairs go into one index per strain.

Outputs land in `results/`: per-TE copy-number table, per-k-mer depths, baselines, and
per-TE depth track plots.

## Validation

- Unit tests and a synthetic-read test with known copy numbers (`tests/`).
- `scripts/subsample.sh` makes a small read subset for plumbing checks (too shallow for a real baseline).
- Sanity check on real data: mPing should be well above 2x in HEG4 vs NB.

See `docs/ALGORITHM.md` for the full algorithm and its known limits. See `PLAN.md` for the original problem statement. Planned: revisit dividing depth by in-file
multiplicity instead of dropping non-unique k-mers.
