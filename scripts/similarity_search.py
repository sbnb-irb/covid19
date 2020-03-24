# Imports
import os
import pandas as pd
import numpy as np
import collections
import pickle
import h5py
from scipy.spatial.distance import cdist
from scipy.stats import rankdata
from chemicalchecker import ChemicalChecker
from chemicalchecker.core.signature_data import DataSignature
from chemicalchecker.util.parser import Converter

if not os.path.exists("outputs"):
    os.mkdir("outputs")

# Read downloaded data from literature
df = pd.read_csv("data/Literature_Candidates.tsv", delimiter="\t")

# Load CC signatures
cc = ChemicalChecker("/aloy/web_checker/package_cc/paper/")
s3 = cc.signature("ZZ.004", "sign3")
universe = set(s3.keys)
universe_conn = dict([(k.split("-")[0], k) for k in s3.keys])

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
        d[inchikey] += [(r[0], r[3], r[4], r[6], r[7], r[8])]
    except:
        continue
R = []
for k, v in d.items():
    score = np.max([x[1] for x in v])
    descriptions = ";".join([x[2] for x in v if type(x[2]) is str])
    links1 = [x[3] for x in v if type(x[3]) is str]
    links2 = [x[4] for x in v if type(x[4]) is str]
    links3 = [x[5] for x in v if type(x[5]) is str]
    links  = links1 + links2 + links3
    links  = ";".join(links)
    R += [(k, k, score, descriptions, links)]
df_lit = pd.DataFrame(R, columns = ["inchikey", "name", "score", "descriptions", "links"])
fn = "outputs/df_lit.pkl"
with open(fn, "wb") as f:
    pickle.dump(df_lit, f)

# Get vectors
iks_lit, V_lit = s3.get_vectors(list(df_lit["inchikey"]))
scores = {}
for r in df_lit[["inchikey", "score"]].values:
    scores[r[0]] = r[1]
scs_lit = np.array([scores[k] for k in iks_lit])
iks_can, V_can = s3.keys, s3[:]
scs_can = []
for k in iks_can:
    if k in scores:
        scs_can += [scores[k]]
    else:
        scs_can += [-2]
scs_can = np.array(scs_can)

# Metadata
isdrug_can = np.array([1]*len(iks_can)) # provisional
nam_can = iks_can # provisional
nam_lit = iks_lit # provisional

# Compute similarities
dist_euc = cdist(V_can, V_lit, metric="euclidean")
dist_cos = cdist(V_can, V_lit, metric="cosine")
def tautize(V):
    W = np.zeros(V.shape)
    for j in range(0, V.shape[1]):
        ranks  = rankdata(-V[:,j], method="ordinal")
        W[:,j] = ranks / np.max(ranks)
    return W
tau_euc = tautize(dist_euc)
tau_cos = tautize(dist_cos)
tau = (tau_euc + tau_cos) / 2
with h5py.File("outputs/dist.h5", "w") as hf:
    hf.create_dataset("euc" , data = dist_euc)
    hf.create_dataset("cos" , data = dist_cos)
    hf.create_dataset("tau" , data = tau)
    hf.create_dataset("rows", data = np.array(iks_can, DataSignature.string_dtype()))
    hf.create_dataset("cols", data = np.array(iks_lit, DataSignature.string_dtype()))
    hf.create_dataset("scs_rows", data = scs_can)
    hf.create_dataset("scs_cols", data = scs_lit)
    hf.create_dataset("isdrug_rows", data = isdrug_can)

# Get matching scores for the candidate molecules
def similarities_as_matrix(tau, min_score):
    if min_score is not None:
        mask = scs_lit >= min_score
        tau  = tau[:,mask]
    sim_avg = np.mean(tau, axis=1)
    sim_q66 = np.percentile(tau, 66, axis=1)
    cou_95  = np.sum(tau > 0.95, axis=1)
    cou_99  = np.sum(tau > 0.99, axis=1)
    cou_999 = np.sum(tau > 0.999, axis=1)
    max_idx = np.argmax(tau, axis = 1)
    M = np.vstack([sim_avg, sim_q66, cou_95, cou_99, cou_999, max_idx]).T
    return M

def similarities(tau, min_score, sort_by="count_99"):
    if min_score is not None:
        if min_score < 0 or min_score > 3:
            raise Exception("min_score must be None, 0, 1, 2 or 3")
    M = similarities_as_matrix(tau, min_score)
    results = {
        "inchikey": iks_can,
        "name": nam_can,
        "is_drug": isdrug_can,
        "lit_score": scs_can,
        "sim_avg": M[:,0],
        "sim_q66": M[:,1],
        "count_95": M[:,2],
        "count_99": M[:,3],
        "count_999": M[:,4],
        "best_inchikey": [iks_lit[int(idx)] for idx in M[:,5]],
        "best_name": [nam_lit[int(idx)] for idx in M[:,5]]
    }
    df = pd.DataFrame(results)
    if min_score is None:
        fn = "df_cand_all.pkl"
    else:
        fn = "df_cand_sc%d.pkl" % min_score
    fn = os.path.join("outputs", fn)
    with open(fn, "wb") as f:
        pickle.dump(df, f)
    return df

similarities(tau, min_score=None)
similarities(tau, min_score=0)
similarities(tau, min_score=1)
similarities(tau, min_score=2)
similarities(tau, min_score=3)

print("Done!")