import os
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import numpy as np
import collections
import h5py
import pandas as pd
from matplotlib import cm
import matplotlib as mpl
import matplotlib.patches as patches
from chemicalchecker.util.plot.diagnosticsplot import set_style, coord_color
from sklearn.neighbors import NearestNeighbors
from scipy.stats import gaussian_kde
from MulticoreTSNE import MulticoreTSNE as TSNE
from tqdm import tqdm
from scipy.stats import fisher_exact
set_style()

script_path = os.path.dirname(os.path.realpath(__file__))
output_path = os.path.join(script_path, '..', 'web', 'static', 'images', 'docu')

OUTPATH = output_path

def literature_figure(simtype):

    SIMTYPE = simtype

    df_lit = pd.read_csv(os.path.join(script_path, "../web/data/df_lit_%s.csv" % SIMTYPE), sep="\t")
   
    evids = ["Failed in clinics", "Text mining", "Computational", "Preclinical", "Clinics", "Clinics COVID19"]
    moas  = ["Unknown", "Host factor", "Virus entry", "Protease inh.", "RNA trans./rep.", "Immunomodulator"]

    eabb = ["Fail", "TM", "Comp", "Pre", "Clin", "CoV"]
    mabb = ["UNK", "HF", "VE", "PI", "RNA", "IM"]

    cap = df_lit.shape[0]/3

    def lit_barplot(ax, df, column):
        counts = collections.defaultdict(int)
        x = []
        y = []
        if column == "Evidence":
            cats = evids
        else:
            cats = moas
        z = []
        for i, cat in enumerate(cats):
            x += [i]
            y += [len(df[df[column] == cat])]
            z += [cat]
        norm = mpl.colors.Normalize(vmin=0,vmax=cap)
        cmap = cm.get_cmap("Spectral")
        color = cmap(norm(y))
        if column == "Evidence":
            for i in range(0, len(x)):
                rect = patches.Rectangle((x[i]-0.2, 0), 0.5, y[i], facecolor=color[i], edgecolor="black", lw=1)
                ax.add_patch(rect)
            ax.set_xticks(x)
            ax.set_xticklabels(eabb)
            ax.set_ylabel("Counts")
            ax.set_ylim(0, np.max(y)*1.05)
            ax.set_xlim(-0.5, len(x)-0.5)
            ax.set_title("Evidence")
        else:
            for i in range(0, len(x)):
                rect = patches.Rectangle((0, x[i]-0.2), y[i], 0.5, facecolor=color[i], edgecolor="black", lw=1)
                ax.add_patch(rect)
            ax.set_yticks(x)
            ax.set_yticklabels(mabb)
            ax.set_xlabel("Counts")
            ax.set_xlim(0, np.max(y)*1.05)
            ax.set_ylim(-0.5, len(x)-0.5)
            ax.set_title("MoA")
        
    def lit_overlaps(ax, df):
        x = []
        y = []
        z = []
        c = []
        for i, e in enumerate(evids):
            for j, m in enumerate(moas):
                a  = len(df[(df["Evidence"] == e) & (df["MoA"] == m)])
                b  = len(df[(df["Evidence"] == e) | (df["MoA"] == m)])
                z += [a/b]
                c += [a]
                x += [i]
                y += [j]
        x = np.array(x)
        y = np.array(y)
        norm = mpl.colors.Normalize(vmin=0,vmax=cap)
        cmap = cm.get_cmap("Spectral")
        c = cmap(norm(c))
        z = np.array(z)
        ax.scatter(x,y,color=c,s=z*2000, edgecolor="black", lw=0.8)
        ax.set_xticks(sorted(set(x)))
        ax.set_yticks(sorted(set(y)))
        ax.set_xticklabels(eabb)
        ax.set_yticklabels(mabb)
        ax.set_xlim(-0.5, len(set(x))-0.5)
        ax.set_ylim(-0.5, len(set(y))-0.5)
        ax.set_title("Overlap")

    def lit_nn(ax):
        with h5py.File(os.path.join(script_path, "../web/data/dist_%s.h5" % SIMTYPE), "r") as hf:
            V = hf["V_lit_trim"][:]
        nn = NearestNeighbors(6)
        nn.fit(V)
        dist = nn.kneighbors(V)[0][:,1:]
        ks = [1,2,3,4,5]
        norm = mpl.colors.Normalize(1,5)
        cmap = cm.get_cmap("Spectral")
        color = cmap(norm(ks))
        for i,k in enumerate(ks):
            scores = dist[:,i]
            kernel = gaussian_kde(scores)
            x = np.linspace(0, np.max(scores), 1000)
            y = kernel(x)
            ax.plot(x,y,color=color[i],label="k=%d" % (k))
        ax.set_title("NN intra literature")
        ax.set_xlabel("Distance")
        ax.set_ylabel("Density")
        ylim = ax.get_ylim()
        ax.set_ylim(0,ylim[1])
        ax.legend()

    with h5py.File(os.path.join(script_path, "../web/data/dist_%s.h5" % SIMTYPE), "r") as hf:
        print(hf.keys())
        evid = hf["evi_cols"][:]
        ranks = hf["ranks"][:]
        ranks_w = hf["ranks_w"][:]
        support = hf["support"][:]

    evi_colordict = {
        0: coord_color("E"),
        1: coord_color("A"),
        2: coord_color("B"),
        3: coord_color("C"),
        4: coord_color("D")
    }
    moa_colordict = {
        0: "gray",
        1: coord_color("E"),
        2: coord_color("A"),
        3: coord_color("B"),
        4: coord_color("C"),
        5: coord_color("D")        
    }

    def image(ax, ranks, evid, max_n):
        idxs = np.argsort(-support)
        M = ranks[idxs][:max_n]
        idxs = np.argsort(-evid)
        M = M[:,idxs]
        evid = evid[idxs]
        evid[evid < 0] = 0
        ax.set_xlabel("Literature drugs")
        ax.set_ylabel("Candidates (ranked by support)")
        x = np.array([i for i in range(0, M.shape[1])])
        y = np.array([j for j in range(0, M.shape[0])])
        norm = mpl.colors.Normalize(vmin=0,vmax=3)
        cmaps = {4: cm.get_cmap("Greens"),
                 3: cm.get_cmap("Blues"),
                 2: cm.get_cmap("Purples"),
                 1: cm.get_cmap("Reds"),
                 0: cm.get_cmap("Oranges")}
        for j in tqdm(range(0, M.shape[1])):
            r = M[:,j]
            ys = y[r > 0]
            cmap = cmaps[evid[j]]
            color = cmap(norm(r[r>0]))
            ax.scatter([j]*len(ys), ys, color=color, s=0.2)
        ax.axvline(np.sum(evid>=4)-0.5, color=evi_colordict[4], lw=1, label="CoV")
        ax.axvline(np.sum(evid>=3)-0.5, color=evi_colordict[3], lw=1, label="Clin")
        ax.axvline(np.sum(evid>=2)-0.5, color=evi_colordict[2], lw=1, label="Pre")
        ax.axvline(np.sum(evid>=1)-0.5, color=evi_colordict[1], lw=1, label="Comp")
        ax.axvline(np.sum(evid>=0)+0.5, color=evi_colordict[0], lw=1, label="TM")
        ax.set_ylim(M.shape[0], 0)
        ax.set_xlim(0, M.shape[1])
        ax.set_title("Similarity matrix")
        ax.legend()

    fig = plt.figure(constrained_layout=True, figsize=(8.5,5))
    gs = fig.add_gridspec(2, 4)
    ax = fig.add_subplot(gs[0,0])
    lit_barplot(ax, df_lit, column="Evidence")
    ax = fig.add_subplot(gs[1,1])
    lit_barplot(ax, df_lit, column="MoA")
    ax = fig.add_subplot(gs[0,1])
    lit_nn(ax)
    ax = fig.add_subplot(gs[1,0])
    lit_overlaps(ax, df_lit)
    ax = fig.add_subplot(gs[2:])
    image(ax, ranks, evid, 2500)
    plt.savefig(os.path.join(OUTPATH, "literature_%s.png" % SIMTYPE), dpi=300)


def projections_figure(simtype):
    SIMTYPE = simtype

    with h5py.File(os.path.join(script_path, "../web/data/dist_%s.h5" % SIMTYPE), "r") as hf:
        V_can_trim = hf["V_can_trim"][:]
        iks_can_trim = hf["iks_can_trim"][:]
        V_lit      = hf["V_lit_trim"][:]
        iks_lit    = hf["cols"][:]
        evi_lit    = hf["evi_cols"][:]
        moa_lit    = hf["moa_cols"][:]
        V = np.vstack([V_lit, V_can_trim])
        
    mani = TSNE()
    Pcan = mani.fit_transform(V)

    with h5py.File(os.path.join(script_path, "../web/data/dist_%s.h5" % SIMTYPE), "r") as hf:
        V_rnd_trim = hf["V_rnd_trim"][:]
        iks_can_trim = hf["iks_rnd_trim"][:]
        V = np.vstack([V_lit, V_rnd_trim])
        
    mani = TSNE()
    Prnd = mani.fit_transform(V)

    def projection_moa(ax, P, title, legend):
        P_lit = P[:len(moa_lit)]
        P_can = P[len(moa_lit):]
        ax.scatter(P_can[:,0], P_can[:,1], color="lightgray", s=1)
        moa_colors = {0: "gray",
                      1: coord_color("A"),
                      2: coord_color("B"),
                      3: coord_color("C"),
                      4: coord_color("D"),
                      5: coord_color("E")}
        moa_legend = {0: "UNK",
                      1: "HF",
                      2: "VF",
                      3: "PI",
                      4: "RNA",
                      5: "IM"}
        for c in [0, 1, 2, 3, 4, 5]:
            mask = moa_lit == c
            ax.scatter(P_lit[mask,0], P_lit[mask,1], color=moa_colors[c], edgecolor="white", label=moa_legend[c])
        ax.set_xlabel("TSNE1")
        ax.set_ylabel("TSNE2")
        ax.set_title(title)
        if legend:
            ax.legend()
        
    def projection_evi(ax, P, title, legend):
        evi_lit[evi_lit < 0] = 0
        P_lit = P[:len(evi_lit)]
        P_can = P[len(evi_lit):]
        ax.scatter(P_can[:,0], P_can[:,1], color="lightgray", s=1)
        evi_colors = {0: coord_color("A"),
                      1: coord_color("B"),
                      2: coord_color("C"),
                      3: coord_color("D")}
        evi_legend = {0: "Comp",
                      1: "Preclin",
                      2: "Clin",
                      3: "COV19"}
        for c in [0, 1, 2, 3]:
            mask = evi_lit == c
            ax.scatter(P_lit[mask,0], P_lit[mask,1], color=evi_colors[c], edgecolor="white", label=evi_legend[c])
        ax.set_xlabel("TSNE1")
        ax.set_ylabel("TSNE2")
        ax.set_title(title)
        if legend:
            ax.legend()    

    fig = plt.figure(constrained_layout=True, figsize=(6,6))
    gs = fig.add_gridspec(2,2)
    ax = fig.add_subplot(gs[0,0])
    projection_moa(ax, Prnd, title="MoA | 10k CC random", legend=True)
    ax = fig.add_subplot(gs[0,1])
    projection_evi(ax, Prnd, title="Evidence | 10k CC random", legend=True)
    ax = fig.add_subplot(gs[1,0])
    projection_moa(ax, Pcan, title="MoA | 10k cand. (by support)", legend=False)
    ax = fig.add_subplot(gs[1,1])
    projection_evi(ax, Pcan, title="Evidence | 10k cand. (by support)", legend=False)
    plt.savefig(os.path.join(OUTPATH, "projections_%s.png" % SIMTYPE), dpi=300)


def atc_figure(simtype):

    def atc_test(ax, simtype=simtype, atc_level="B", top=100, exclude_self=True, yaxis=True):
        with h5py.File(os.path.join(script_path, "../web/data/dist_%s.h5" % simtype), "r") as hf:
            iks_lit = hf["cols"][:]
            iks_can = hf["rows"][:]
            support = hf["support"][:]
        idxs = np.argsort(-support)
        support = support[idxs]
        iks_can = iks_can[idxs]
        pos_iks = iks_can[:top]
        neg_iks = iks_can[top:]
        def get_atcs(level="B"):
            atcs = collections.defaultdict(set)
            with open(os.path.join(script_path, "../data/atcs.tsv"), "r") as f:
                for l in f:
                    l = l.rstrip("\n").split("\t")
                    atcs[l[0].split("-")[0]].update([x.split(":")[1] for x in l[1].split(",") if "%s:" % atc_level in x])
            return atcs
        atcs = get_atcs()    
        atc_descr = {
            "A": "Alim. tract & Met.",
            "B": "Blood etc.",
            "C": "Cardiovascular",
            "D": "Dermatol.",
            "G": "Genitourinary etc.",
            "H": "Syst. horm",
            "J": "Antiinfectives",
            "L": "Antineoplastic",
            "M": "Musculoskeletal",
            "N": "Nervous",
            "P": "Antiparasitic",
            "R": "Respiratory",
            "S": "Sensory organs",
            "V": "Various"
        }
        descr_sort = sorted([x for x in atc_descr.keys()])
        universe = set(atcs.keys())
        if exclude_self:
            universe.difference(set([x.split("-")[0] for x in iks_lit]))
        Pos = set()
        Neg = set()
        for ik in pos_iks:
            ik = ik.split("-")[0]
            if ik not in universe:
                continue
            Pos.update([ik])
        for ik in neg_iks:
            ik = ik.split("-")[0]
            if ik not in universe:
                continue
            Neg.update([ik])
        Neg = Neg.difference(Pos)
        universe = universe.intersection(Pos.union(Neg))
        atc_sets = collections.defaultdict(set)
        for k,v in atcs.items():
            if k not in universe:
                continue
            for x in v:
                atc_sets[x].update([k])
        atc_cats = sorted(atc_sets.keys())
        x = []
        y = []
        enriched = []
        for atc in atc_cats:
            v = atc_sets[atc]
            A = len(v.intersection(Pos))
            B = len(Pos) - A
            C = len(v.intersection(Neg))
            D = len(universe) - (A+B+C)
            odds, pval = fisher_exact([[A,B],[C,D]], alternative="greater")
            x += [odds]
            y += [-np.log10(pval)]
        x = np.array(x)
        y = np.array(y)
        labs = np.array(atc_cats)
        x_ = [i for i in range(0, len(y))]
        if top == 100:
            color = coord_color("A")
        if top == 1000:
            color = coord_color("E")
        if top == 10000:
            color = coord_color("D")
        ax.scatter(y, x_, s=np.sqrt(x*1000)+1, color=color)
        for i in range(0, len(y)):
            ax.plot([0, y[i]], [x_[i], x_[i]], color=color)
        ax.set_yticks(x_)
        if yaxis:
            ax.set_yticklabels([atc_descr[l] for l in labs])
        else:
            ax.set_yticklabels(labs)
        ax.set_xlabel("-log10 P-value")
        ax.set_title("ATC in top %d" % top)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[1], ylim[0])
        
    fig, axs = plt.subplots(1,3, figsize=(6, 3.5))
    axs = axs.flatten()
    atc_test(axs[0], simtype=simtype, atc_level="A", top=100, yaxis=True)
    atc_test(axs[1], simtype=simtype, atc_level="A", top=1000, yaxis=False)
    atc_test(axs[2], simtype=simtype, atc_level="A", top=10000, yaxis=False)
    plt.tight_layout()
    plt.savefig(os.path.join(script_path, "../web/static/images/docu/atc_%s.png" % simtype), dpi=300)


def significance_figure(simtype):
    pass


def candidates(simtype):

    def load_df(simtype, evidence, moa):
        if evidence is None or evidence < 0:
            evi_suf = "eviall"
        else:
            evi_suf = "evi%d" % evidence
        if moa is None:
            moa_suf = "moaall"
        else:
            moa_suf = "moa%d" % moa
        fn = "../web/data/df_cand_%s_%s_%s.csv" % (simtype, evi_suf, moa_suf)
        df = pd.read_csv(fn, sep="\t")
        return df

    def candidate_view(df):
        return df[["inchikey", "name", "is_drug", "support", "lpv_5", "lpv_4", "lpv_3", "top1_name", "top2_name", "top3_name"]][df["evidence"]==-2]

    
def do_plots(simtype):
    literature_figure(simtype)
    projections_figure(simtype)
    atc_figure(simtype)
    significance_figure(simtype)
    candidates(simtype)

if __name__ == "__main__":
    do_plots("cc")
    do_plots("fp")

