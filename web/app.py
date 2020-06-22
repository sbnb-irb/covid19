import os
import pandas as pd
from flask import Flask, request, render_template, jsonify, send_file

from .tables import BaseDataTables, ServerSideTable

app = Flask(__name__)
app.debug = True

app_path = os.path.dirname(os.path.realpath(__file__))


def get_literature_data():
    file_path = os.path.join(app_path, 'data', 'df_lit_cc.csv')
    data = pd.read_csv(file_path, sep="\t")
    data.fillna('!N/A', inplace=True)

    def urls_to_html(links):
        links = links.split(',')
        html = ''
        for idx, link in enumerate(links, 1):
            html += '<a target="_blank" href=%s>[%i] </a>' % (link[1:-1], idx)
        return html

    # link is the column with hyperlinks
    data['References'] = data['References'].apply(urls_to_html)
    print('LOADED %s' % file_path)
    print('LENGTH %s' % len(data))
    print('COLUMNS %s' % str(data.columns))
    print('EXAMPLE\n%s' % str(data.iloc[0]))
    return data


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
    data.fillna('!N/A', inplace=True)
    print('LOADED %s' % file_path)
    print('LENGTH %s' % len(data))
    print('COLUMNS %s' % str(data.columns))
    print('EXAMPLE\n%s' % str(data.iloc[0]))
    return data


def get_query_data(query_id, signature):
    file_path = os.path.join(app_path, 'static', 'images', 'docu',
                             'query_%s_%s.csv' % (query_id, signature))
    data = pd.read_csv(file_path, sep="\t")
    data.fillna('!N/A', inplace=True)
    print('LOADED %s' % file_path)
    print('LENGTH %s' % len(data))
    print('COLUMNS %s' % str(data.columns))
    print('EXAMPLE\n%s' % str(data.iloc[0]))
    return data


df_candidates = get_candidate_data()
df_candidates_collection = df_candidates.to_dict(orient='records')


def get_candidate_update(signature='cc', evidence='', moa=''):
    filen_name = 'similarity_update.txt'
    file_path = os.path.join(app_path, 'data', filen_name)
    with open(file_path, 'r') as fh:
        last_update = fh.readline()
    return last_update


@app.route('/')
def index():
    return render_template('index.html', columns=df_candidates.columns,
                           last_update=get_candidate_update())


@app.route('/literature')
def literature():
    df_literature = get_literature_data()
    return render_template('literature.html', columns=df_literature.columns)


@app.route('/docs')
def docs():
    df = get_query_data('1', 'cc')
    return render_template('docs.html', columns=df.columns,
                           last_update=get_candidate_update())


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/_literature_table')
def get_literature_table():
    df_literature = get_literature_data()
    df_literature_collection = df_literature.to_dict(orient='records')
    columns = df_literature.columns
    collection = df_literature_collection
    results = BaseDataTables(request, columns, collection).output_result()
    print("RESULTS", str(results))
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
    print("RESULTS", str(results))
    return jsonify(results)


@app.route('/_table_title')
def get_table_title():
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


@app.route('/_table_query')
def get_query_table():
    print('get_query_table', request.args)
    query_id = request.args.get('query')
    signature = request.args.get('signature')
    df_query = get_query_data(query_id, signature)
    df_query_collection = df_query.to_dict(orient='records')
    columns = list()
    for order, col in enumerate(df_query.columns, 1):
        col_dict = {
            "data_name": col,
            "column_name": col,
            "default": "",
            "order": order,
            "searchable": True
        }
        columns.append(col_dict)
    collection = df_query_collection
    results = ServerSideTable(request, columns, collection).output_result()
    print("RESULTS", str(results))
    return jsonify(results)


@app.route('/download')
def download_file():
    print('download_file', request.args)
    signature = request.args.get('signature')
    evidence = request.args.get('evidence')
    moa = request.args.get('moa')
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
    return send_file(file_path, as_attachment=True)


if __name__ == '__main__':
    app.run()
