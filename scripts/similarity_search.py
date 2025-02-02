# Imports
import os
import h5py
import shutil
import gspread
import datetime
import collections
import numpy as np
import pandas as pd
import csv
from glob import glob
from tqdm import tqdm
from oauth2client.service_account import ServiceAccountCredentials

from chemicalchecker import ChemicalChecker
from chemicalchecker.util.parser import Converter
from chemicalchecker.core.signature_data import DataSignature

from plots import do_plots

# TODO make this an argparse
script_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(script_path, '..', 'data')
similarity_path = os.path.join(script_path, '..', 'tmp', 'webdata')
plot_path = os.path.join(script_path, '..', 'tmp', 'docu')
final_similarity_path = os.path.join(script_path, '..', 'web', 'data')
final_plot_path = os.path.join(
    script_path, '..', 'web', 'static', 'images', 'docu')

os.makedirs(similarity_path, exist_ok=True)
os.makedirs(plot_path, exist_ok=True)

MAX_ROWS = 10000


def get_raw_literature():
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(
        os.path.join(data_path, 'covid19_repo.json'), scope)
    client = gspread.authorize(creds)
    # get all worksheets
    sheets = client.open('Drugs against SARS-CoV-2').worksheets()
    # check if candidates sheet is available
    for sheet in sheets:
        if sheet.title == 'Candidates':
            break
    if sheet.title != 'Candidates':
        titles = [s.title for s in sheets]
        raise Exception('"Candidates" sheet not found! Available: %s' % titles)
    # get row as pandas dataframe
    records = sheet.get_all_records()
    df = pd.DataFrame(records)
    # check if expected column names are there
    if not all(df.columns[:3] == ['Fazy', 'Drug Name', 'SMILES']):
        raise Exception('Wrong headers! Found: %s' % str(list(df.columns)))
    '''
    updated = sheet.updated
    print('LITERATURE UPDATED:', updated)
    dest_file = os.path.join(similarity_path, "literature_update.txt")
    with open(dest_file, 'wb') as fh:
        fh.write(updated)
    '''
    print(df)
    return df


evidence_legend = {
    -2: "NA",
    -1: "Failed in clinics",
    0: "Text mining",
    1: "Computational",
    2: "Preclinical",
    3: "Clinics",
    4: "Clinics COVID19"
}
moa_legend = {
    -2: "NA",
    -1: "Not given",
    0: "Unknown",
    1: "Host factor",
    2: "Virus entry",
    3: "Protease inh.",
    4: "RNA trans./rep.",
    5: "Immunomodulator"
}


def main(simtype):

    # Load CC signatures
    cc = ChemicalChecker("/aloy/web_checker/package_cc/paper/")
    if simtype == "cc":
        print("Loading CC signatures")
        neig = cc.signature("ZZ.004", "neig3")
        sign = cc.signature("ZZ.004", "sign3")
    else:
        print("Loading Morgan fingerprints")
        neig = cc.signature("A1.001", "neig1")
        sign = cc.signature("A1.001", "sign1")
    print("Getting universe")
    universe = set(neig.row_keys)
    universe_conn = dict([(k.split("-")[0], k) for k in neig.row_keys])

    print('Fetching literature annotation.')
    literature_file = os.path.join(similarity_path, 'literature.csv') 
    try:
        df = get_raw_literature()
        # merge references to single column
        df['References'] = df[[
            'Reference', 'Reference (2)',
            'Reference (3)', 'Reference (4)']].agg('","'.join, axis=1)
        df['References'] = '"' + df['References'] + '"'
        df['References'] = df['References'].str.replace(',""', '')
        df.to_csv(literature_file, index=False, sep="\t")
    except Exception as ex:
        print('Exception while fetching data:')
        print(str(ex))
        print('Using cached version')
        df = pd.read_csv(literature_file, sep="\t")

    print('Reading text-mining candidates')
    df_tm = pd.read_csv(os.path.join(data_path, "textmining_selection.tsv"),
                        header=None, names=["inchikey"], delimiter="\t")

    # Drugbank
    with open(os.path.join(data_path, "db2ikey.tsv"), "r") as f:
        druglist = []
        druglist_conn = []
        favconn_names = {}
        favconn_iks = {}
        for l in f:
            l = l.rstrip("\n").split("\t")
            druglist += [l[-1]]
            druglist_conn += [l[-1].split("-")[0]]
            favconn_names[l[-1].split("-")[0]] = l[1]
            favconn_iks[l[-1].split("-")[0]] = l[-1]
        druglist = set(druglist)
        druglist_conn = set(druglist_conn)

    # Names
    ik_name = {}
    with open(os.path.join(data_path, "inchikeys_names.csv"), "r") as f:
        for l in f:
            l = l.rstrip("\n")
            ik = l[:27]
            name = l[28:]
            ik_name[ik] = name

    # Smiles
    ik_smiles = {}
    with open(os.path.join(data_path, "inchikeys_smiles.tsv"), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            ik_smiles[r[0]] = r[1]

    print("Converting inchikeys of literature")
    # Convert to inchikeys
    conv = Converter()
    d = collections.defaultdict(list)
    for r in df.values:
        try:
            mol = conv.smiles_to_inchi(r[2])
            inchikey, inchi = mol
            smi = conv.inchi_to_smiles(inchi)
            if inchikey not in universe:
                if inchikey.split("-")[0] not in universe_conn:
                    continue
                else:
                    inchikey = universe_conn[inchikey.split("-")[0]]
            d[inchikey] += [(r[0], r[1], r[3], r[5], r[4], r[12])]
            ik_smiles[inchikey] = smi
        except Exception as ex:
            print(ex)
            continue
    print("Adding textmining")
    for i, ik in enumerate(list(df_tm["inchikey"])):
        if ik not in universe:
            if ik.split("-")[0] not in universe_conn:
                continue
            else:
                ik = universe_conn[ik.split("-")[0]]
        d[ik] += [("tm%d" % i, ik, 0, -1, "Text mining", "", "", "")]
    print("Assembling")
    R = []
    for k, v in d.items():
        c = k.split("-")[0]
        if c in favconn_names:
            name = favconn_names[c]
        else:
            if k in ik_name:
                name = ik_name[k]
            else:
                name = k
        evi = [x[2] for x in v if not np.isnan(x[2])]
        if len(evi) > 0:
            evi = np.max(evi)
        else:
            evi = -2
        moa = [x[3] for x in v if not np.isnan(x[3])]
        if len(moa) > 0:
            moa = int(np.max(moa))
        else:
            moa = -2
        urls = v[0][5]
        R += [(k, name, evi, moa, urls)]
    print("Aggregate literature candidates")
    d_lit = collections.defaultdict(list)
    for r in R:
        d_lit[r[0].split("-")[0]] += [r]
    R = []
    for k, v in d_lit.items():
        if len(v) == 1:
            R += [v[0]]
        else:
            evids_ = [x[2] for x in v]
            idx = np.argmax(evids_)
            n = v[idx][1]
            e = v[idx][2]
            m = v[idx][3]
            u = v[idx][4]
            if k in favconn_iks:
                ik = favconn_iks[k]
            else:
                ik = v[idx][0]
            R += [(ik, n, e, m, u)]
    R_ = []
    for r in R:
        R_ += [(r[0], r[1], r[2], evidence_legend[r[2]], moa_legend[r[3]], r[4])]
    df_lit = pd.DataFrame(
        R_, columns=["InChIKey", "Name", "Level", "Evidence", "MoA", "References"])
    df_lit = df_lit.sort_values(["Level", "InChIKey"], ascending=[False, True])
    print("...saving literature for web")
    dest_file = os.path.join(similarity_path, "df_lit_%s.csv" % simtype)
    df_lit.to_csv(dest_file, index=False, sep="\t")
    print("...keeping tractable version of literature")
    df_lit = pd.DataFrame(
        R, columns=["inchikey", "name", "evidence", "moa", "references"])
    print(df_lit.shape)
    print("Getting precomputed similarities")

    # Get precomputed similarities
    iks_can = neig.row_keys
    iks_can_set = set(iks_can)
    iks_lit = np.array(
        [k for k in sorted(df_lit["inchikey"]) if k in iks_can_set])
    iks_can_idxs = dict((k, i) for i, k in enumerate(iks_can))
    idxs_lit = [iks_can_idxs[k] for k in iks_lit]
    with h5py.File(neig.data_path, "r") as hf:
        nn = hf["indices"][idxs_lit]

    print("Keeping evidence")
    evids = {}
    for r in df_lit[["inchikey", "evidence"]].values:
        evids[r[0]] = r[1]
    evi_lit = np.array([evids[k] for k in iks_lit])
    evi_can = []
    for k in iks_can:
        if k in evids:
            evi_can += [evids[k]]
        else:
            evi_can += [-2]
    evi_can = np.array(evi_can)

    print("Keeping moa")
    moas = {}
    for r in df_lit[["inchikey", "moa"]].values:
        moas[r[0]] = r[1]
    moa_lit = np.array([moas[k] for k in iks_lit])
    moa_can = []
    for k in iks_can:
        if k in moas:
            moa_can += [moas[k]]
        else:
            moa_can += [-2]
    moa_can = np.array(moa_can)

    # Metadata
    print("Metadata: is drug")
    # Drugs
    isdrug_can = []
    for k in iks_can:
        if k.split("-")[0] in druglist_conn:
            isdrug_can += [1]
        else:
            isdrug_can += [0]
    isdrug_can = np.array(isdrug_can)

    print("Metadata: names")
    # Names
    nam_can = []
    for k in iks_can:
        if k.split("-")[0] in favconn_names:
            nam_can += [favconn_names[k.split("-")[0]]]
        else:
            if k in ik_name:
                nam_can += [ik_name[k]]
            else:
                nam_can += [k]
    nam_lit = []
    for k in iks_lit:
        if k.split("-")[0] in favconn_names:
            nam_lit += [favconn_names[k.split("-")[0]]]
        else:
            if k in ik_name:
                nam_lit += [ik_name[k]]
            else:
                nam_lit += [k]

    print("Reorganizing ranks")
    ranks_raw = np.full((len(iks_can), len(iks_lit)), -1, dtype=np.int32)
    vals = [i for i in range(0, nn.shape[1])]
    for j in range(0, ranks_raw.shape[1]):
        ranks_raw[nn[j, :], j] = vals

    # clean
    del nn

    print("...applying p-value cutoffs (1e-5=3, 1e-4=2, 1e-3=1)")
    ranks = np.array(ranks_raw, dtype=np.int32)
    cutoffs = {
        "lpv_5": int(np.round(ranks.shape[0] * 1e-5, 0)),
        "lpv_4": int(np.round(ranks.shape[0] * 1e-4, 0)),
        "lpv_3": int(np.round(ranks.shape[0] * 1e-3, 0))
    }
    print(cutoffs)
    ranks[ranks == -1] = 100000
    ranks[np.logical_and(ranks < cutoffs["lpv_5"], ranks >= 0)] = -3
    ranks[np.logical_and(ranks < cutoffs["lpv_4"], ranks >= 0)] = -2
    ranks[np.logical_and(ranks < cutoffs["lpv_3"], ranks >= 0)] = -1
    ranks[ranks >= 0] = 0
    ranks = ranks * (-1)
    print(ranks.shape)
    print("...weighting ranks by evidence (<=0=1, 1=2, 2=3, 3=4, 4=5)")
    evis = np.array(evi_lit, dtype=np.int8)
    evis[evis < 0] = 0
    evis = evis + 1
    ranks_w = ranks * evis
    support = np.sum(ranks_w, axis=1)
    print("...aggregate by inchikey connectivity")
    iks_can_conn_idx = collections.defaultdict(list)
    for i, ik in enumerate(iks_can):
        iks_can_conn_idx[ik.split("-")[0]] += [i]
    keep_idx_ = []
    iks_can_ = []
    evi_can_ = []
    moa_can_ = []
    nam_can_ = []
    isdrug_can_ = []
    conn_lit_idx = {}
    assert (len(iks_lit) == len(idxs_lit))
    for i, k in enumerate(iks_lit):
        conn_lit_idx[k.split("-")[0]] = idxs_lit[i]
    for conn, idxs in iks_can_conn_idx.items():
        if len(idxs) == 1:
            idx = idxs[0]
            keep_idx_ += [idx]
            iks_can_ += [iks_can[idx]]
            evi_can_ += [evi_can[idx]]
            moa_can_ += [moa_can[idx]]
            nam_can_ += [nam_can[idx]]
            isdrug_can_ += [isdrug_can[idx]]
        else:
            if conn in conn_lit_idx:
                idx = conn_lit_idx[conn]
            else:
                supps = [support[idx] for idx in idxs]
                idx = idxs[np.argmax(supps)]
            # keep inchikey
            keep_idx_ += [idx]
            # inchikey
            if conn in favconn_iks:
                iks_can_ += [favconn_iks[conn]]
            else:
                iks_can_ += [iks_can[idx]]
            # highest evidence
            evi_can_ += [np.max([evi_can[i] for i in idxs])]
            # highest moa (arbitrary)
            moa_can_ += [np.max([moa_can[i] for i in idxs])]
            # name
            if conn in favconn_names:
                nam_can_ += [favconn_names[conn]]
            else:
                nam_can_ += [nam_can[idx]]
            # conservative is drug
            isdrug_can_ += [np.max([isdrug_can[i] for i in idxs])]
    sort_idxs = np.argsort(iks_can_)
    iks_can = np.array(iks_can_)[sort_idxs]
    keep_idx = np.array(keep_idx_, dtype=np.int32)[sort_idxs]
    evi_can = np.array(evi_can_, dtype=np.int8)[sort_idxs]
    moa_can = np.array(moa_can_, dtype=np.int8)[sort_idxs]
    nam_can = [nam_can_[i] for i in sort_idxs]
    isdrug_can = np.array(isdrug_can_, dtype=np.int8)[sort_idxs]
    ranks = ranks[keep_idx]
    ranks_raw = ranks_raw[keep_idx]
    ranks_w = ranks_w[keep_idx]
    support = support[keep_idx]
    assert (ranks_w.shape[0] == ranks_raw.shape[0] == ranks.shape[0] ==
            len(isdrug_can) == len(evi_can) == len(iks_can) == len(nam_can) ==
            len(moa_can))
    del iks_can_
    del keep_idx
    del keep_idx_
    del evi_can_
    del moa_can_
    del nam_can_
    del isdrug_can_

    print("Getting vectors for plotting")
    iks_lit_trim, V_lit_trim = sign.get_vectors(iks_lit)
    iks_can_trim, V_can_trim = sign.get_vectors(
        iks_can[np.argsort(-support)][:10000])
    rnd_idxs = np.arange(len(iks_can))
    np.random.shuffle(rnd_idxs)
    iks_rnd_trim, V_rnd_trim = sign.get_vectors(
        iks_can[rnd_idxs][:10000])
    print("...saving")
    dest_file = os.path.join(similarity_path, "dist_%s.h5" % simtype)
    with h5py.File(dest_file, "w") as hf:
        hf.create_dataset("ranks_raw", data=ranks_raw)
        hf.create_dataset("ranks", data=ranks)
        hf.create_dataset("ranks_w", data=ranks_w)
        hf.create_dataset("support", data=support)
        hf.create_dataset("rows", data=np.array(
            iks_can, DataSignature.string_dtype()))
        hf.create_dataset("cols", data=np.array(
            iks_lit, DataSignature.string_dtype()))
        hf.create_dataset("evi_rows", data=evi_can)
        hf.create_dataset("evi_cols", data=evi_lit)
        hf.create_dataset("moa_rows", data=moa_can)
        hf.create_dataset("moa_cols", data=moa_lit)
        hf.create_dataset("isdrug_rows", data=isdrug_can)
        hf.create_dataset("iks_lit_trim", data=np.array(
            iks_lit_trim, DataSignature.string_dtype()))
        hf.create_dataset("V_lit_trim", data=V_lit_trim)
        hf.create_dataset("iks_can_trim", data=np.array(
            iks_can_trim, DataSignature.string_dtype()))
        hf.create_dataset("V_can_trim", data=V_can_trim)
        hf.create_dataset("iks_rnd_trim", data=np.array(
            iks_rnd_trim, DataSignature.string_dtype()))
        hf.create_dataset("V_rnd_trim", data=V_rnd_trim)
        hf.create_dataset("nam_cols", data=np.array(
            nam_lit, DataSignature.string_dtype()))

    # Get matching scores for the candidate molecules
    def scoring(ranks, min_evidence, moa):
        if min_evidence is None:
            if moa is None:
                mask = None
            else:
                mask = moa_lit == moa
        else:
            mask1 = evi_lit >= min_evidence
            if moa is None:
                mask = mask1
            else:
                mask2 = moa_lit == moa
                mask = np.logical_and(mask1, mask2)
        if not np.any(mask):
            mask = None
        if mask is None:
            mask = np.array([True for _ in range(0, ranks.shape[1])])
        support = np.sum(ranks_w[:, mask], axis=1, dtype=np.int32)
        keep = np.argsort(-support)[:MAX_ROWS]
        keep = keep[support[keep] > 0]
        support = support[keep]
        cou_pv5 = np.sum(ranks[keep][:, mask] >= 3, axis=1, dtype=np.int16)
        cou_pv4 = np.sum(ranks[keep][:, mask] >= 2, axis=1, dtype=np.int16)
        cou_pv3 = np.sum(ranks[keep][:, mask] >= 1, axis=1, dtype=np.int16)
        top1_idx = []
        top2_idx = []
        top3_idx = []
        mask_idxs = np.argwhere(mask).ravel()
        for idx in tqdm(keep):
            worth_idxs = np.argwhere(ranks[idx, :] > 0).ravel()
            worth_idxs = np.intersect1d(worth_idxs, mask_idxs)
            vals = ranks_raw[idx, worth_idxs]
            top_idxs = worth_idxs[np.argsort(vals)[:3]]
            n = len(top_idxs)
            if n == 1:
                top1_idx += [top_idxs[0]]
                top2_idx += [-1]
                top3_idx += [-1]
            elif n == 2:
                top1_idx += [top_idxs[0]]
                top2_idx += [top_idxs[1]]
                top3_idx += [-1]
            else:
                top1_idx += [top_idxs[0]]
                top2_idx += [top_idxs[1]]
                top3_idx += [top_idxs[2]]
        top1_idx = np.array(top1_idx, dtype=np.int32)
        top2_idx = np.array(top2_idx, dtype=np.int32)
        top3_idx = np.array(top3_idx, dtype=np.int32)
        M = np.vstack([support, cou_pv5, cou_pv4, cou_pv3,
                       top1_idx, top2_idx, top3_idx]).astype(np.int32)
        M = M.T
        return M, keep, ranks.shape[0], np.sum(mask)

    def similarities(simtype, min_evidence, moa):
        if min_evidence is not None:
            if min_evidence < 0 or min_evidence > 4:
                raise Exception("min_evidence must be None, 0, 1, 2 or 4")
        M, keep, n_rows, n_cols = scoring(ranks, min_evidence, moa)
        iks_lit_ = list(iks_lit) + [""]
        nam_lit_ = list(nam_lit) + [""]
        results = {
            "inchikey": iks_can[keep],
            "name": [nam_can[i] for i in keep],
            "is_drug": isdrug_can[keep],
            "evidence": evi_can[keep],
            "moa": moa_can[keep],
            "support": M[:, 0],
            "lpv_5": M[:, 1],
            "lpv_4": M[:, 2],
            "lpv_3": M[:, 3],
            "top1_inchikey": [iks_lit_[int(idx)] for idx in M[:, 4]],
            "top2_inchikey": [iks_lit_[int(idx)] for idx in M[:, 5]],
            "top3_inchikey": [iks_lit_[int(idx)] for idx in M[:, 6]],
            "top1_name": [nam_lit_[int(idx)] for idx in M[:, 4]],
            "top2_name": [nam_lit_[int(idx)] for idx in M[:, 5]],
            "top3_name": [nam_lit_[int(idx)] for idx in M[:, 6]]
        }
        df = pd.DataFrame(results)
        print("...adding smiles")
        smiles = []
        for ik in list(df["inchikey"]):
            if ik in ik_smiles:
                smiles += [ik_smiles[ik]]
            else:
                smiles += [np.nan]
        df["smiles"] = smiles
        print("...filtering")
        if simtype == "cc":
            caption = "CC similarities"
        else:
            caption = "Chemical similarities"
        if min_evidence is None:
            evi_suf = "eviall"
        else:
            evi_suf = "evi%d" % min_evidence
        if moa is None:
            moa_suf = "moaall"
        else:
            moa_suf = "moa%d" % moa
        caption += " against %d drugs from the COVID19 literature" % n_cols
        fn = "df_cand_%s_%s_%s.csv" % (simtype, evi_suf, moa_suf)
        df.to_csv(os.path.join(similarity_path, fn), index=False, sep="\t")
        return fn, caption

    print("Saving legend")
    legend = list()
    for min_evidence in [None, 0, 1, 2, 3, 4]:
        for moa in [None, 0, 1, 2, 3, 4, 5]:
            print("Sim type %s | Min evidence %s | MoA %s" %
                  (simtype, str(min_evidence), str(moa)))
            fn, caption = similarities(
                simtype=simtype, min_evidence=min_evidence, moa=moa)
            legend.append({
                'evidence': min_evidence,
                'moa': moa,
                'filename': fn,
                'caption': caption})
    legend = pd.DataFrame(legend)
    dest_file = os.path.join(similarity_path, "legend_%s.csv" % simtype)
    legend.to_csv(dest_file, index=False, sep="\t")

    print("Doing plots")
    try:
        do_plots(simtype, data_path, similarity_path, plot_path)
    except Exception:
        print("Plotting failed!")


if __name__ == "__main__":
    main(simtype="cc")
    main(simtype="fp")
    print("Copying SIMILARITY files to final destination.")
    for src in glob(os.path.join(similarity_path, '*.csv')):
        print('%s to %s' % (src, final_similarity_path))
        shutil.copy(src, final_similarity_path)
    for src in glob(os.path.join(similarity_path, "dist_*.h5")):
        print('%s to %s' % (src, final_similarity_path))
        shutil.copy(src, final_similarity_path)
    print("Copying PLOT files to final destination.")
    for src in glob(os.path.join(plot_path, '*.png')):
        print('%s to %s' % (src, final_plot_path))
        shutil.copy(src, final_plot_path)
    for src in glob(os.path.join(plot_path, '*.csv')):
        print('%s to %s' % (src, final_plot_path))
        shutil.copy(src, final_plot_path)

    now = datetime.datetime.now()
    updated = str(now.strftime("%c"))
    print('UPDATE DONE:', updated)
    dest_file = os.path.join(final_similarity_path, "similarity_update.txt")
    with open(dest_file, 'w') as fh:
        fh.write(updated)
