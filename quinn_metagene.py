#!/usr/bin/env python3

#### module load python/3.9.19

#### USAGE: script.py co.bed dsb.bed 500 50 100

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================================
# usage system & input files
# =====================================================

GENE_FILE = "/home2/qyj2/from_maizegdb/v5/b73_v5_gene_tss_tts_sort.bed"

CO_FILE  = sys.argv[1]
DSB_FILE = sys.argv[2]

# independent bin sizes
FLANK_BIN  = int(sys.argv[3]) 
BODY_BIN   = int(sys.argv[4]) 
INNER_BIN  = int(sys.argv[5])

# region sizes
FLANK_SIZE = 10000
INNER_SIZE = 2000

SMOOTH = 5

# =====================================================
# define gene 
# =====================================================

n_up    = FLANK_SIZE // FLANK_BIN
n_inner = INNER_SIZE // INNER_BIN
n_body  = BODY_BIN  

TOTAL = n_up + n_inner + n_body + n_inner + n_up

# offsets
OFF_TSS  = n_up
OFF_BODY = OFF_TSS + n_inner
OFF_TTS  = OFF_BODY + n_body
OFF_DOWN = OFF_TTS + n_inner

# =====================================================
# load & process data
# =====================================================

genes = pd.read_csv(
    GENE_FILE, sep="\t", header=None,
    names=["chr","start","end","len","dot","strand","type"]
)

co = pd.read_csv(
    CO_FILE, sep="\t", header=None,
    names=["chr","start","end","count"]
)
dsb = pd.read_csv(
    DSB_FILE, sep="\t", header=None,
    names=["chr","start","end","count"]
)

co["mid"]  = ((co.start + co.end)//2).astype(int)
dsb["mid"] = ((dsb.start + dsb.end)//2).astype(int)

co_chr  = {c:df for c,df in co.groupby("chr")}
dsb_chr = {c:df for c,df in dsb.groupby("chr")}

# =====================================================
# gene profile function
# =====================================================

def gene_profile(g, signal_dict):

    profile = np.zeros(TOTAL)
    if g.chr not in signal_dict:
        return profile

    sig = signal_dict[g.chr]
    m = sig.mid.values
    c = sig["count"].values

    # strand normalization
    if g.strand == "+":
        TSS, TTS = g.start, g.end
    else:
        TSS, TTS = g.end, g.start

    # gene body limits
    body_start = TSS + INNER_SIZE
    body_end   = TTS - INNER_SIZE
    body_len   = max(1, body_end - body_start)

    # --------------- upstream ---------------
    mask = (m >= TSS-FLANK_SIZE) & (m < TSS)
    idx = ((m[mask]-(TSS-FLANK_SIZE))//FLANK_BIN).astype(int)
    valid = (idx >= 0) & (idx < n_up)
    np.add.at(profile, idx[valid], c[mask][valid])

    # --------------- TSS inner ---------------
    mask = (m >= TSS) & (m < TSS+INNER_SIZE)
    idx = OFF_TSS + ((m[mask]-TSS)//INNER_BIN).astype(int)
    valid = (idx >= OFF_TSS) & (idx < OFF_BODY)
    np.add.at(profile, idx[valid], c[mask][valid])

    # --------------- scaled gene body ---------------
    mask = (m >= body_start) & (m < body_end)
    scaled = ((m[mask]-body_start)/body_len) * n_body
    idx = OFF_BODY + scaled.astype(int)
    valid = (idx >= OFF_BODY) & (idx < OFF_TTS)
    np.add.at(profile, idx[valid], c[mask][valid])

    # --------------- TTS inner ---------------
    mask = (m >= TTS-INNER_SIZE) & (m < TTS)
    idx = OFF_TTS + ((m[mask]-(TTS-INNER_SIZE))//INNER_BIN).astype(int)
    valid = (idx >= OFF_TTS) & (idx < OFF_DOWN)
    np.add.at(profile, idx[valid], c[mask][valid])

    # --------------- downstream ---------------
    mask = (m >= TTS) & (m < TTS+FLANK_SIZE)
    idx = OFF_DOWN + ((m[mask]-TTS)//FLANK_BIN).astype(int)
    valid = (idx >= OFF_DOWN) & (idx < TOTAL)
    np.add.at(profile, idx[valid], c[mask][valid])

    return profile

# =====================================================
# aggregation function
# =====================================================

def aggregate(signal_dict):

    profiles = []

    for _, g in genes.iterrows():
        p = gene_profile(g, signal_dict)
        if p.sum() == 0:
            continue
        profiles.append(p)

    profiles = np.vstack(profiles)
    return profiles.mean(axis=0) #metagene profile

# =====================================================
# smoothing function
# =====================================================

def smooth(x, w=3):
    return np.convolve(x, np.ones(w)/w, mode="same")

# =====================================================
# run aggregation & smoothing
# =====================================================

co_raw  = aggregate(co_chr)
dsb_raw = aggregate(dsb_chr)

co_raw  = smooth(co_raw, SMOOTH)
dsb_raw = smooth(dsb_raw, SMOOTH)

# -----------------------------------------------------
# plot: define boundaries
# -----------------------------------------------------

TSS_START   = OFF_TSS
TSS_END     = OFF_BODY

BODY_START  = OFF_BODY
BODY_END    = OFF_TTS

TTS_START   = OFF_TTS
TTS_END     = OFF_DOWN

BODY_CENTER = (BODY_START + BODY_END) / 2

# -----------------------------------------------------
# plot: helper function to format x-axis
# -----------------------------------------------------

def style_dual_anchor(ax):

    # Shade TSS and TTS windows
    ax.axvspan(TSS_START, TSS_END, alpha=0.08)
    ax.axvspan(TTS_START, TTS_END, alpha=0.08)

    # Anchor TSS and TTS
    ax.axvline(TSS_START, linestyle="--", linewidth=1)
    ax.axvline(TTS_END, linestyle="--", linewidth=1)

    # Clean x-axis labels
    xticks = [
        0,
        TSS_START,
        BODY_CENTER,
        TTS_END,
        TOTAL - 1
    ]

    xlabels = [
        "-10 kb",
        "TSS",
        "Gene body",
        "TTS",
        "+10 kb"
    ]

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_xlim(0, TOTAL - 1)


# -----------------------------------------------------
# plot: generate metagene
# -----------------------------------------------------
fig, ax = plt.subplots(1, 1, figsize=(6, 4))

x = np.arange(TOTAL)

ax.plot(x, co_raw, label="CO")
ax.plot(x, dsb_raw, label="DSB")
ax.set_title("Dual-anchor metagene")
ax.set_ylabel("Density")
ax.legend()

style_dual_anchor(ax)

plt.tight_layout()
plt.savefig("dual_anchor_metagene.pdf", bbox_inches="tight")
plt.show()




