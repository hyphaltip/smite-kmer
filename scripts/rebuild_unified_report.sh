#!/usr/bin/env bash
# Rebuild the unified per-strain tables (k=21 and k=31) and the report in docs/report.
# Strains: the earlier panel, the landraces and the stratified 3K list, whichever have a histogram and a depth table.
# usage (sbatch): sbatch scripts/rebuild_unified_report.sh
set -euo pipefail
root=${PROJ_ROOT:?set PROJ_ROOT}
source /etc/profile.d/modules.sh
module load jellyfish/2.3.0
cd "$root"
python3 - <<'PY'
import os
import pandas as pd
base = ["HEG4","NB","B031","IRIS_313-11418","CX55","IRIS_313-10761","CX269","IRIS_313-10871","IRIS_313-8321",
        "IRIS_313-8385","IRIS_313-11652","IRIS_313-8171","A119","A123","A157","NB_DRR028132"]
stream = pd.read_csv("resources/Chen2019/stream100.tsv", sep="\t").strain.tolist()
for k in (21, 31):
    cand = list(dict.fromkeys(base + stream))
    ok = [s for s in cand if os.path.exists(f"work/{s}.k{k}.histo") and os.path.exists(f"work/{s}.k{k}.depth.tsv")]
    pd.DataFrame({"strain": ok, "reads": "work"}).to_csv(f"work/unified_k{k}.tsv", sep="\t", index=False)
    print(f"k{k}: {len(ok)} strains of {len(cand)} candidates")
PY
for k in 21 31; do
  mkdir -p results/unified/k$k
  python3 scripts/mitekmer.py analyze --fasta RiceTE.fa --samples work/unified_k$k.tsv --workdir work --k $k \
    --reference NB --out results/unified/k$k --skip-stats
done
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/pandoc/2.17.1.1/bin:$PATH
python3 scripts/make_te_panel_report.py --k21 results/unified/k21 --k31 results/unified/k31 --out docs/report \
  --tex-bin /opt/linux/rocky/8.x/x86_64/pkgs/texlive/20220403/bin/x86_64-linux
