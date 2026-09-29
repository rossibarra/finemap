"""Plot cM/Mb from Ogut v5 and finemap_hierarchical_v5 over a 10 Mb region, with a gene rug.

Usage: python scripts/plot_random_region.py [--seed N | --region Chr1:158600000]
"""
import argparse
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int)
ap.add_argument("--region", help="CHROM:START (10 Mb window starting at START)")
args = ap.parse_args()

hier = pd.read_csv("data/finemap_hierarchical_v5.bed", sep="\t", header=None,
                   names=["chr","start","end","cM0","cM1","rate"])
og = pd.read_csv("data/ogut_v5.csv").sort_values(["chr","pos_v5"])
gff = pd.read_csv("data/v5.genes.gff3", sep="\t", header=None, comment="#",
                  usecols=[0, 2, 3, 4], names=["chr","feature","start","end"])
genes = gff[gff.feature == "gene"]
lens = hier.groupby("chr").end.max()

if args.region:
    chrom, start = args.region.split(":"); start = int(start); label = args.region
else:
    seed = args.seed if args.seed is not None else int(np.random.default_rng().integers(1e6))
    rng = np.random.default_rng(seed)
    chrom = rng.choice(lens.index, p=lens/lens.sum())
    start = int(rng.integers(0, (lens[chrom]-10_000_000)//100_000))*100_000
    label = f"random seed {seed}"
end = start + 10_000_000

h = hier[(hier.chr==chrom)&(hier.end>start)&(hier.start<end)]
o = og[og.chr==chrom]
# Ogut cM interpolated at the same 100 kb edges, differenced -> cM/Mb
edges = np.arange(start, end+1, 100_000)
cm = np.interp(edges, o.pos_v5, o.cM)
orate = np.diff(cm)/0.1
mk = o.pos_v5[(o.pos_v5>=start)&(o.pos_v5<end)]
g = genes[(genes.chr.str.lower()==chrom.lower())&(genes.end>start)&(genes.start<end)]
gmid = (g.start + g.end)/2

fig, (ax, rug) = plt.subplots(2, 1, figsize=(11,5), sharex=True,
                              gridspec_kw={"height_ratios":[10,1], "hspace":0.05})
ax.stairs(orate, edges/1e6, label=f"Ogut v5 (linear interp., 100 kb; {len(mk)} markers)", color="#d95f02", lw=1.6)
ax.stairs(h.rate.values, np.r_[h.start.values, h.end.values[-1]]/1e6, label="finemap hierarchical (100 kb)", color="#1b9e77", lw=1.6)
for p in mk: ax.axvline(p/1e6, ymin=0, ymax=0.025, color="#d95f02", lw=0.6)
ax.set_ylabel("cM/Mb"); ax.set_ylim(bottom=0); ax.legend(frameon=False)
ax.set_title(f"{chrom}:{start/1e6:.1f}-{end/1e6:.1f} Mb  ({label})")
for s in ["top","right"]: ax.spines[s].set_visible(False)

rug.vlines(gmid/1e6, 0, 1, color="0.2", lw=0.5, alpha=0.5)
rug.set_ylim(0, 1); rug.set_yticks([]); rug.set_ylabel(f"genes\n(n={len(g)})", rotation=0, ha="right", va="center")
for s in ["top","right","left"]: rug.spines[s].set_visible(False)
rug.set_xlim(start/1e6, end/1e6); rug.set_xlabel(f"{chrom} position (Mb, B73 v5)")

fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.11)
out = f"results/region_{chrom}_{start/1e6:.1f}-{end/1e6:.1f}Mb_ogut_vs_hier.png"
fig.savefig(out, dpi=150); print(out)
print("Ogut cM:", cm[-1]-cm[0], " hier cM:", h.cM1.iloc[-1]-h.cM0.iloc[0], " genes:", len(g))
