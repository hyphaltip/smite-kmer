# Handoff: smite-kmer

Written 2026-10-26 for a change of HPC account. Read this file first. It lists what exists, what is running, and what to do next.

## Goal

Test whether k-mer depth in short reads recovers transposable-element (TE) copy numbers. The main test is mPing. The paper is Chen et al. 2019 (Nat Commun 10:641, DOI 10.1038/s41467-019-08451-3). Its Supplementary Data 1 gives read-depth mPing copy numbers for about 3,000 rice accessions.

## Repository

- Path: `/bigdata/stajichlab/jstajich/projects/smite-kmer` (also `/rhome/jstajich/projects/smite-kmer`, the same repo).
- Remote: `https://github.com/hyphaltip/smite-kmer.git`, branch `main`.
- Last pushed commit: `fd305d6`. The working tree is clean, except for files listed under "Uncommitted" below.
- Gitignored and not in git: `reads/`, `work/`, `results/`, `logs/`, `sra_cache/`, `fastq_tmp/`, and `*.jf`.

## Method (summary)

Full detail is in `docs/ALGORITHM.md`.

1. Count canonical k-mers per accession with jellyfish, at k=21 and k=31, quality filter Q20.
2. Take the TE k-mers from `RiceTE.fa`, keeping those that occur once in the file.
3. Single-copy depth (baseline) = highest histogram peak above the error valley. Returns NaN when no clear peak exists (`baseline_peak`).
4. Copy number = median(TE k-mer depth) / baseline.
5. Compare to the reference strain NB (Nipponbare, 52 copies in the paper).

## Key results so far

- Low-mPing panel (9 strains with a baseline), k=21: Pearson r 0.969, Spearman 0.979, median ratio to paper 1.12.
- NB: 52.4 copies (paper 52). CX140 (paper's Nipponbare entry, streamed from ENA): 56.6.
- Success set (76 accessions plotted at k=21, 66 at k=31): Spearman 0.981 (k=21) and 0.978 (k=31). The median ratio to paper is 1.08 to 1.23 per bin up to 100 copies. The highest plotted estimate is B160 (paper 180): 119 at k=21, 125 at k=31.
- No plotted accession reaches the burst range (231 to 503). HEG4 (485), A119 (360) and A123 (276) are in the range, but the paper gives no per-line values for them.
- Reports (dated 2026-10-26):
  - Main: `docs/report/REPORT.md` and `REPORT.pdf`, built from `results/success/k21`.
  - Success set: `docs/report/success/REPORT.md` and `REPORT.pdf`.

## Pipeline scripts

| Script | Use |
|---|---|
| `scripts/count.sh` | Count one accession from local FASTQ folders (folder discovery, symlinks followed). |
| `scripts/profile_stream.sh` | Count one or more k from read files (local, URL, or truncated-gzip prefix). Writes the histogram, depth table, stats and GC-bias table. |
| `scripts/sra_node_task.sh` | Current route. One accession on one compute node: prefetch from NCBI SRA into `$SCRATCH`, convert with parallel-fastq-dump, check every run's pair count against ENA, count k-mers, write tables to `work/`. |
| `scripts/sra_node_array_task.sh` | SLURM array wrapper for `sra_node_task.sh`. |
| `scripts/sra_fetch.sh`, `scripts/sra_pipeline.sh` | Login-node prefetch route. Superseded by `sra_node_task.sh`, which keeps all downloads on the compute node. |
| `scripts/submit.sh`, `scripts/submit_stream.sh` | Older submitters. Defaults are account `stajichlab`, partition `stajichlab`. |
| `scripts/stream_task.sh` | Old ENA streaming array task. Superseded. |
| `scripts/download_panel.sh` | Old aria2c ENA download. Superseded. |
| `scripts/mitekmer.py` | `kmers`, `query`, `baseline`, `gcbias`, `analyze`, `profile` subcommands. |
| `scripts/make_success_figure.py` | Figure and table for the success set. |
| `scripts/make_prelim_report.py` | Main report (Markdown, PDF, figures). Needs `TEX_BIN` for the PDF. |
| `scripts/select_panel.py`, `scripts/select_stratified.py` | Choose accessions from the paper table. |
| `scripts/make_sra_runs.py` | Build a run list (strain, run, ENA read count). |

Rule for read sources: one source per accession. ENA is used for metadata only, and NCBI SRA for reads.

## Data and reference files

- `resources/Chen2019/41467_2019_8451_MOESM4_ESM.xlsx`: paper Supplementary Data 1.
- `resources/Chen2019/PRJEB6180.ena_runs.tsv.gz`: ENA run table, with read counts.
- `resources/Chen2019/stream100.tsv`: 100 stratified accessions (paper mPing bins), with runs.
- `resources/Chen2019/lowmping.panel.tsv`: the 10 low-mPing panel accessions.
- `resources/Chen2019/rerun_candidates.tsv`, `rerun34_strains.tsv`, `sra_runs.rerun35.tsv`: the 35 accessions being rerun, and their runs.
- `resources/strain_groups.tsv`: each strain's group (burst, reference, low panel, landrace) and paper value.
- `samples.tsv`, `samples.all.tsv`, `samples.positive.tsv`, `samples.landraces.tsv`: sample sheets.
- Local reads: `reads/HEG4`, `reads/NB` (symlinks to `/bigdata/wesslerlab/shared/Rice/`), `reads/lowmping/`, `reads/landraces/`, `reads/nipponbare_extra/`.
- The NB_DRR028132 folder holds only the first 2.27 GB of R1, from an incomplete download. Use it only with the caveat in the report.

## Results locations

- `work/<strain>.k<k>.histo`, `.depth.tsv`, `.stats`, `.jf` (the jf index is deleted after streaming runs).
- `work/gcbias/<strain>.k<k>.gcbias.tsv`: finished marker and GC-bias table. An empty table means no baseline.
- `results/success/k21/`, `results/success/k31/`: analysis output (copy numbers, baseline, track plots) for the successful accessions.
- `results/full/`: older 12-strain analysis.
- `work/stream100_ready.tsv`, `work/analyze_k21.tsv`, `work/analyze_k31.tsv`: lists of which accessions have complete outputs.

## Running now

Array job **30818516**, account `stajichlab`, partition `stajichlab`. It runs `sra_node_array_task.sh` over `resources/Chen2019/rerun34_strains.tsv` (34 accessions), at most 4 at a time. Each task uses 8 CPUs, 64 GB and 4 hours. At the last check, one task had completed, four were running, and the rest were pending on the group CPU limit.

Logs: `logs/nodetask.rerun34.30818516_<task>.out`.

The 35th accession, **IRIS_313-8114**, is the test. It finished in job 30811372. Its output is in `work/`. It peaked at 31.6 GB of memory.

Check progress with:

    squeue -j 30818516
    sacct -j 30818516 -X --format=JobID,State,Elapsed,MaxRSS

## Next steps

1. When array 30818516 finishes, check each task's log for `ok, <n> pairs` lines, a `mismatch` line, or an OOM kill.
2. For each accession, confirm these exist in `work/`: `k21.histo`, `k21.depth.tsv`, `k31.histo`, `k31.depth.tsv`, and `gcbias/<strain>.k31.gcbias.tsv`.
3. Rerun `scripts/mitekmer.py analyze` on the complete sample sheet for k=21 and k=31, into a new output directory. Use `--reference NB`. The bootstrap step makes this slow (it ran for more than 40 minutes per k on a login node earlier, and did not finish in the first attempt). Run it as a SLURM job, not on the login node. Runtime was not measured on the cluster.
4. Regenerate `docs/report/success/` with `scripts/make_success_figure.py`, then the report with a new date. Commit the outputs.
5. Regenerate `docs/report/REPORT.md` if needed, with `scripts/make_prelim_report.py --results <analysis dir> --k 21`.

## Open decisions and known issues

- **Storage rule:** The global CLAUDE.md says to keep SRA and FASTQ files off node-local scratch. The user directed the per-node route to use `$SCRATCH`. Downloads and FASTQ are lost when a job ends. Only the count tables are kept.
- **Memory:** A 64 GB request is adequate for the largest accession so far (31.6 GB peak). The earlier stream failures had no logged cause. Some of them match out-of-memory kills, but this was not confirmed for every failure.
- **Failed earlier streams:** The original 100-accession array (30321736) had 47 failed tasks. Of these, 12 already had complete tables. The other 35 are the rerun above.
- **Baseline:** Some accessions have no clear single-copy peak. They are reported as "no baseline" and have no copy number. CX55 and IRIS_313-10761 (k=31) are examples. Copy numbers for them are missing, not zero.
- **GC correction:** Implemented in `mitekmer.py`, but not applied to any reported result.
- **Paper values:** The paper gives no per-line values for A119, A123 and A157, so they are not compared to the paper.
- **Shallow accessions:** Strains at about 5x to 9x depth have the weakest baselines (for example CX55 and IRIS_313-8385).

## Environment notes

- Modules do not persist between shell calls. To use jellyfish, run `module load jellyfish/2.3.0` in the same command. Or set `PATH` to `/opt/linux/rocky/8.x/x86_64/pkgs/jellyfish/2.3.0/bin` and `LD_LIBRARY_PATH` to its `lib` directory.
- For PDFs, set `TEX_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux` and put pandoc on `PATH` (`/opt/linux/rocky/8.x/x86_64/pkgs/pandoc/2.17.1.1/bin`).
- SRA tools: `module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7` (`sra_node_task.sh` does this).
- `sbatch --wrap` runs under `sh`. Put Bash-only logic in a script file.
- The `sbatch` default account and partition in the scripts is now `stajichlab`. Override with `ACCOUNT` or `PARTITION`.

## Uncommitted

- `docs/HANDOFF.md` (this file), until it is committed.

## Contacts

Project owner: Jason Stajich (jstajich).
