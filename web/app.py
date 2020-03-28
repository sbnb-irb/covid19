import os
import pandas as pd
from flask import Flask, request, render_template, jsonify

from .tables import BaseDataTables, ServerSideTable

app = Flask(__name__)
app.debug = True

app_path = os.path.dirname(os.path.realpath(__file__))


def get_literature_data():
    file_path = os.path.join(app_path, 'data', 'df_lit_cc.csv')
    data = pd.read_csv(file_path, sep="\t")
    print('LOADED %s' % file_path)
    print('LENGTH %s' % len(data))
    return data


df_literature = get_literature_data()
df_literature_collection = df_literature.to_dict(orient='records')


def get_candidate_data(signature='cc', evidence='', moa=''):
    if evidence == '':
        evidence = 'all'
    else:
        evidence = int(evidence)
    if moa == '':
        moa = 'all'
    else:
        moa = int(moa)
    filen_name = 'df_cand_%s_evi%s_moa%s.csv' % (signature, evidence, moa)
    file_path = os.path.join(app_path, 'data', filen_name)
    data = pd.read_csv(file_path, sep="\t")
    data['name'] = data['name'].fillna('')
    print('LOADED %s' % file_path)
    print('LENGTH %s' % len(data))
    return data


df_candidates = get_candidate_data()
df_candidates_collection = df_candidates.to_dict(orient='records')


@app.route('/')
def index():
    return render_template('index.html', columns=df_candidates.columns)


@app.route('/literature')
def literature():
    return render_template('literature.html', columns=df_literature.columns)


@app.route('/_literature_table')
def get_literature_table():
    columns = df_literature.columns
    collection = df_literature_collection
    results = BaseDataTables(request, columns, collection).output_result()
    return jsonify(results)


@app.route('/_candidate_table')
def get_candidate_table():
    print('get_candidate_table', request.args)
    signature = request.args.get('signature')
    evidence = request.args.get('evidence')
    moa = request.args.get('moa')
    df_candidates = get_candidate_data(signature, evidence, moa)
    df_candidates_collection = df_candidates.to_dict(orient='records')
    columns = list()
    for order, col in enumerate(df_candidates.columns, 1):
        col_dict = {
            "data_name": col,
            "column_name": col,
            "default": "",
            "order": order,
            "searchable": True
        }
        columns.append(col_dict)
    collection = df_candidates_collection
    results = ServerSideTable(request, columns, collection).output_result()
    return jsonify(results)


@app.route('/_table_title')
def get_table_title(signature='cc', evidence='', moa=''):
    print('get_table_title', request.args)
    signature = request.args.get('signature')
    evidence = request.args.get('evidence')
    moa = request.args.get('moa')
    filen_name = 'legend_%s.csv' % signature
    file_path = os.path.join(app_path, 'data', filen_name)
    df = pd.read_csv(file_path, sep="\t")
    if evidence == '':
        df = df[(df['evidence'].isnull())]
    else:
        evidence = float(evidence)
        df = df[(df['evidence'] == evidence)]
    if moa == '':
        df = df[(df['moa'].isnull())]
    else:
        moa = float(moa)
        df = df[(df['moa'] == moa)]
    print('get_table_title RESULT', df)
    assert(len(df) == 1)
    return jsonify(df.iloc[0]['caption'])


if __name__ == '__main__':
    app.run()
