import numpy as np
import wget
import os
import zipfile
import gzip
import shutil

ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"
path = "/aloy/scratch/mduran/pubchem_bioassays/"

def download():
    indices = np.arange(1, 1349000, 1000)
    for i in range(0, len(indices)-1):
        a = "%d" % indices[i]
        a = a.zfill(7)
        b = "%d" % (indices[i+1]-1)
        b = b.zfill(7)
        fn = "%s_%s.zip" % (a,b)
        fn_inp = os.path.join(ftp_path, fn)
        fn_out = os.path.join(path, fn)
        if os.path.exists(fn_out):
            continue
        print(fn)
        try:
            wget.download(fn_inp, out=fn_out)
        except:
            print(fn, "not found")

def parse_one(fn):
    if os.path.exists(os.path.join(path, fn.split(".")[0]+".tsv")):
        return
    with zipfile.ZipFile(os.path.join(path, fn), 'r') as zip_ref:
        zip_ref.extractall(os.path.join(path, "tmp"))
    outpath = os.path.join(path, "tmp", fn.split(".")[0])
    R = []
    if not os.path.exists(outpath):
        return
    for a in os.listdir(outpath):
        with gzip.open(os.path.join(outpath, a), "rb") as f:
            text = str(f.read())
            aid  = text.split("<PC-ID_id>")[1].split("</PC-ID_id>")[0]
            try:
                descr1 = text.split("<PC-AssayDescription_name>")[1].split("</PC-AssayDescription_name")[0]
                descr2 = ".".join([x.split("</PC-AssayDescription_description_E>")[0] for x in text.split("<PC-AssayDescription_description_E>")[1:]])
                descr  = ".".join([descr1, descr2])
            except:
                continue
            R += [(aid, descr.replace("\n", ".").replace("\t", "."))]
    shutil.rmtree(outpath)
    with open(os.path.join(path, fn.split(".")[0]+".tsv"), "w") as f:
        for r in R:
            f.write("%s\t%s\n" % r)

def parse():
    for fn in os.listdir(path):
        if fn[-4:] != ".zip":
            continue
        print(fn)
        parse_one(fn)

if __name__ == "__main__":
    #download()
    parse()