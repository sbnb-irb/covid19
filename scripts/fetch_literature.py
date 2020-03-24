"""Fetch literature data from shared sheet.

This requires the activation of google APIs (drive and sheet) and generation
of credential file.
ref: [https://towardsdatascience.com/accessing-google-spreadsheet-data-using-python-90a5bc214fd2]
"""
import gspread
import pandas as pd
from oauth2client.service_account import ServiceAccountCredentials


def get_raw_literature():
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(
        'covid19_repo.json', scope)
    client = gspread.authorize(creds)
    sheet = client.open('Drugs against SARS-CoV-2').get_worksheet(1)
    records = sheet.get_all_records()
    df = pd.DataFrame(records)
    return df
