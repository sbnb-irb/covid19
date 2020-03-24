"""Fetch literature data from shared sheet.

This requires the activation of google APIs (drive and sheet) and generation
of credential file.
ref: [https://towardsdatascience.com/accessing-google-spreadsheet-data-using-python-90a5bc214fd2]
"""
import os
import gspread
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials

script_path = os.path.dirname(os.path.realpath(__file__))


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


def save_literature(destination):
    df = get_raw_literature()
    print('Saving table to: %s' % destination)
    df.to_csv(destination, index=False)


def main():
    destination = os.path.join(
        script_path, '..', 'web', 'data', 'literature.csv')
    save_literature(destination)


if __name__ == "__main__":
    main()
