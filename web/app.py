import os
import json
import pandas as pd
from flask import Flask, request, render_template

app = Flask(__name__)
app.debug = True

app_path = os.path.dirname(os.path.realpath(__file__))


def get_literature_table():
    return pd.read_csv(os.path.join(app_path, 'data', 'literature.csv'))


class BaseDataTables:

    def __init__(self, request, columns, collection):
        self.columns = columns
        self.collection = collection
        # values specified by the datatable for filtering, sorting, paging
        self.request_values = request.values
        # results from the db
        self.result_data = None
        # total in the table after filtering
        self.cardinality_filtered = 0
        # total in the table unfiltered
        self.cadinality = 0
        self.run_queries()

    def output_result(self):
        output = {}
        # output['sEcho'] = str(int(self.request_values['sEcho']))
        # output['iTotalRecords'] = str(self.cardinality)
        # output['iTotalDisplayRecords'] = str(self.cardinality_filtered)
        aaData_rows = []
        for row in self.result_data:
            aaData_row = []
            for i in range(len(self.columns)):
                print row, self.columns, self.columns[i]
                aaData_row.append(
                    str(row[self.columns[i]]).replace('"', '\\"'))
            aaData_rows.append(aaData_row)
        output['aaData'] = aaData_rows
        return output

    def run_queries(self):
        self.result_data = self.collection
        self.cardinality_filtered = len(self.result_data)
        self.cardinality = len(self.result_data)


@app.route('/')
def index():
    df = get_literature_table()
    return render_template('index.html', columns=df.columns)
    return 'Hello World!'


@app.route('/_server_data')
def get_server_data():

    df = get_literature_table()
    df['Score'].fillna(-1, inplace=True)

    def make_href(val):
        return '<a href="{}">link</a>'.format(val, val)

    df['Reference'] = df['Reference'].apply(make_href)
    columns = df.columns
    collection = [r.to_dict() for _, r in df.iterrows()]

    results = BaseDataTables(request, columns, collection).output_result()
    # return the results as a string for the datatable
    return json.dumps(results)


if __name__ == '__main__':
    app.run()
