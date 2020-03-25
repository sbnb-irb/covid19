"""Fetch literature data from shared sheet.

This requires the activation of google APIs (drive and sheet) and generation
of credential file.
ref: [https://towardsdatascience.com/accessing-google-spreadsheet-data-using-python-90a5bc214fd2]
"""
import os
import gspread
import numpy as np
import collections
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials

from chemicalchecker import ChemicalChecker
from chemicalchecker.util.parser import Converter

script_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(script_path, '..', 'web', 'data')


def get_raw_literature():
    print('Fetching literature annotation.')
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(
        os.path.join(script_path, 'covid19_repo.json'), scope)
    client = gspread.authorize(creds)
    sheet = client.open('Drugs against SARS-CoV-2').get_worksheet(1)
    records = sheet.get_all_records()
    df = pd.DataFrame(records)
    return df


def get_clean_literature():
    df = get_raw_literature()
    print('Cleaning literature annotation.')
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
        except Exception as ex:
            print(ex)
            continue
    R = []
    for k, v in d.items():
        score = np.max([x[1] for x in v])
        descriptions = ";".join([x[2] for x in v if type(x[2]) is str])
        links1 = [x[3] for x in v if type(x[3]) is str]
        links2 = [x[4] for x in v if type(x[4]) is str]
        links3 = [x[5] for x in v if type(x[5]) is str]
        links = links1 + links2 + links3
        links = ";".join(links)
        R += [(k, k, score, descriptions, links)]
    df_lit = pd.DataFrame(
        R, columns=["inchikey", "name", "score", "descriptions", "links"])
    return df_lit


def save_literature(destination):
    df = get_clean_literature()
    print('Saving table to: %s' % destination)
    df.to_csv(destination, index=False)


def main():
    destination = os.path.join(data_path, 'literature.csv')
    save_literature(destination)


if __name__ == "__main__":
    main()
