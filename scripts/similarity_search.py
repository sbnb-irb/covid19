# Imports
import os
import h5py
import gspread
import collections
import numpy as np
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials

from chemicalchecker import ChemicalChecker
from chemicalchecker.util.parser import Converter
from chemicalchecker.core.signature_data import DataSignature


script_path = os.path.dirname(os.path.realpath(__file__))
output_path = os.path.join(script_path, '..', 'web', 'data')
input_path = os.path.join(script_path, '..', 'data')


def get_raw_literature():
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(
        os.path.join(script_path, 'covid19_repo.json'), scope)
    client = gspread.authorize(creds)
    sheet = client.open('Drugs against SARS-CoV-2').get_worksheet(1)
    records = sheet.get_all_records()
    df = pd.DataFrame(records)
    return df


def main(simtype):

    # Load CC signatures
    cc = ChemicalChecker("/aloy/web_checker/package_cc/paper/")
    if simtype == "cc":
        print("Loading CC signatures")
        neig = cc.signature("ZZ.004", "neig3")
    else:
        print("Loading Morgan fingerprints")
        neig = cc.signature("A1.001", "neig1")

    universe = set(neig.row_keys)
    universe_conn = dict([(k.split("-")[0], k) for k in neig.row_keys])

    print('Fetching literature annotation.')
    df = get_raw_literature()
    literature_file = os.path.join(input_path, 'literature.csv')
    df.to_csv(literature_file, index=False)
    literature_file = os.path.join(output_path, 'literature.csv')
    df.to_csv(literature_file, index=False)

    print("Converting inchikeys")
    # Convert to inchikeys
    conv = Converter()
    d = collections.defaultdict(list)
    for r in df.values:
        try:
            mol = conv.smiles_to_inchi(r[2])
            inchikey, inchi = mol
            if inchikey not in universe:
                if inchikey.split("-")[0] not in universe_conn:
                    continue
                else:
                    inchikey = universe_conn[inchikey.split("-")[0]]
            d[inchikey] += [(r[0], r[1], r[3], r[5], r[4], r[8], r[9], r[10])]
        except Exception as ex:
            print(ex)
            continue
    R = []
    for k, v in d.items():
        name = ";".join([x[1] for x in v if type(x[1]) is str])
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
        descriptions = ";".join([x[4] for x in v if type(x[4]) is str])
        links1 = [x[5] for x in v if type(x[5]) is str]
        links2 = [x[6] for x in v if type(x[6]) is str]
        links3 = [x[7] for x in v if type(x[7]) is str]
        links = links1 + links2 + links3
        links = ";".join(links)
        R += [(k, name, evi, moa, descriptions, links)]
    df_lit = pd.DataFrame(
        R, columns=["inchikey", "name", "evidence",
                    "moa", "descriptions", "links"])
    dest_file = os.path.join(output_path, "df_lit_%s.csv" % simtype)
    df_lit.to_csv(dest_file, index=False)

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
    with open("data/druglist.tsv", "r") as f:
        druglist = []
        for l in f:
            druglist += [l.rstrip("\n")]
        druglist_conn = [d.split("-")[0] for d in druglist]
        druglist = set(druglist)
        druglist_conn = set(druglist_conn)
    isdrug_can = []
    for k in iks_can:
        if k.split("-")[0] in druglist_conn:
            isdrug_can += [1]
        else:
            isdrug_can += [0]
    isdrug_can = np.array(isdrug_can)

    print("Metadata: names")
    #  Names
    ik_name = {}
    with open("data/inchikeys_names.csv", "r") as f:
        for l in f:
            l = l.rstrip("\n")
            ik = l[:27]
            name = l[28:]
            ik_name[ik] = name
    nam_can = []
    for k in iks_can:
        if k in ik_name:
            nam_can += [ik_name[k]]
        else:
            nam_can += [k]
    nam_lit = []
    for k in iks_lit:
        if k in ik_name:
            nam_lit += [ik_name[k]]
        else:
            nam_lit += [k]

    print("Reorganizing ranks")
    ranks = np.full((len(iks_can), len(iks_lit)), -1).astype(np.int)
    vals = [i for i in range(0, nn.shape[1])]
    for j in range(0, ranks.shape[1]):
        ranks[nn[j, :], j] = vals

    # Compute similarities
    dest_file = os.path.join(output_path, "dist_%s.h5" % simtype)
    with h5py.File(dest_file, "w") as hf:
        hf.create_dataset("ranks", data=ranks)
        hf.create_dataset("rows", data=np.array(
            iks_can, DataSignature.string_dtype()))
        hf.create_dataset("cols", data=np.array(
            iks_lit, DataSignature.string_dtype()))
        hf.create_dataset("evi_rows", data=evi_can)
        hf.create_dataset("evi_cols", data=evi_lit)
        hf.create_dataset("isdrug_rows", data=isdrug_can)

    # Get matching scores for the candidate molecules
    def similarities_as_matrix(ranks, min_evidence, moa):
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
        if mask is not None:
            ranks = ranks[:, mask]
        ranks[ranks == -1] = 99999
        cou_tp5 = np.sum(ranks <= 4, axis=1).astype(np.int)
        cou_pv5 = np.sum(ranks <= ranks.shape[0] * 1e-5, axis=1).astype(np.int)
        cou_pv4 = np.sum(ranks <= ranks.shape[0] * 1e-4, axis=1).astype(np.int)
        cou_pv3 = np.sum(ranks <= ranks.shape[0] * 1e-3, axis=1).astype(np.int)
        cou_pv2 = np.sum(ranks <= ranks.shape[0] * 1e-2, axis=1).astype(np.int)
        best_idx = np.argmin(ranks, axis=1).astype(np.int)
        M = np.vstack([cou_tp5, cou_pv5, cou_pv4, cou_pv3,
                       cou_pv2, best_idx]).astype(np.int)
        M = M.T
        return M, ranks.shape[0], ranks.shape[1]

    def similarities(simtype, min_evidence, moa, sort_by="lpv_4"):
        if min_evidence is not None:
            if min_evidence < 0 or min_evidence > 3:
                raise Exception("min_evidence must be None, 0, 1, 2 or 3")
        M, n_rows, n_cols = similarities_as_matrix(ranks, min_evidence, moa)
        results = {
            "inchikey": iks_can,
            "name": nam_can,
            "is_drug": isdrug_can,
            "evidence": evi_can,
            "moa": moa_can,
            "top_5": M[:, 0],
            "lpv_5": M[:, 1],
            "lpv_4": M[:, 2],
            "lpv_3": M[:, 3],
            "lpv_2": M[:, 4],
            "best_inchikey": [iks_lit[int(idx)] for idx in M[:, 5]],
            "best_name": [nam_lit[int(idx)] for idx in M[:, 5]]
        }
        df = pd.DataFrame(results)
        print("...filtering")
        iks = set()
        for col in ["top_5", "lpv_5", "lpv_4", "lpv_3", "lpv_2"]:
            df = df.sort_values(col, ascending=False)
            iks.update(list(df[df[col] > 0]["inchikey"][:10000]))
        df = df[df["inchikey"].isin(iks)]
        caption = "%s similarities" % simtype.upper()
        if min_evidence is None:
            evi_suf = "eviall"
        else:
            evi_suf = "evi%d" % min_evidence
        if moa is None:
            moa_suf = "moaall"
        else:
            moa_suf = "moa%d" % moa
        caption += " (%d cand. x %d lit.)" % (n_rows, n_cols)
        fn = "df_cand_%s_%s_%s.csv" % (simtype, evi_suf, moa_suf)
        df = df.sort_values(sort_by, ascending=False)
        df.to_csv(os.path.join(output_path, fn), index=False)
        return fn, caption

    print("Saving legend")
    legend = list()
    for min_evidence in [None, 0, 1, 2, 3]:
        for moa in [None, 0, 1, 2, 3, 4, 5]:
            print("Sim type %s | Min evidence %s | MoA %s" %
                  (simtype, str(min_evidence), str(moa)))
            fn, caption = similarities(
                simtype=simtype, min_evidence=min_evidence, moa=moa)
            legend.append(
                {'evidence': min_evidence, 'moa': moa, 'filename': fn, 'caption': caption})
    legend = pd.DataFrame(legend)
    dest_file = os.path.join(output_path, "legend_%s.csv" % simtype)
    legend.to_csv(dest_file, index=False)


if __name__ == "__main__":
    main(simtype="cc")
    main(simtype="fp")

print("Done!")
