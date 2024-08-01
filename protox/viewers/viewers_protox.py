import os
import webbrowser
import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
from protox.protocols.protocol_protox import ProtChemProtox
import pandas as pd

class ProtChemProtoxViewer(pwviewer.ProtocolViewer):
    """ Viewer for ProtChemProtox Protocol """
    _label = 'View Protox 3.0 Analysis Results'
    _targets = [ProtChemProtox]

    def __init__(self, **args):
        super().__init__(**args)


    def _defineParams(self, form):
        form.addSection(label='View Scores')
        group_json = form.addGroup('CSV Data')
        group_json.addParam('displayCsv', params.LabelParam,
                            label='Open CSV file: ',
                            help='Show CSV data.')


    def _getVisualizeDict(self):
        return {
            'displayCsv': self._render
        }
    
    def _render(self, paramName=None):
        csv_path = self.protocol._getExtraPath("results.csv")
        self._render_results_table(csv_path)

    def _render_results_table(self, csv_path):
        df = self.read_csv_to_dataframe(csv_path)
        if df is not None:
            html_table = self.dataframe_to_html(df)
            self.save_html_file(html_table)
        else:
           print("No data available to display.")

    def read_csv_to_dataframe(self, csv_path):
        try:
            df = pd.read_csv(csv_path, sep='\t', usecols=lambda col: col != 'Unnamed: 0')
            print(df)
            return df
        except pd.errors.ParserError as e:
            print(f"Error reading CSV: {e}")
            return None

    def dataframe_to_html(self, df):
        html_table = df.to_html(classes='table table-striped', index=False)
        styled_html = f"""
        <style>
            .table {{
                width: 100%;
                border-collapse: collapse;
            }}
            .table th, .table td {{
                text-align: center;
                padding: 8px;
                border: 1px solid #ddd;
            }}
        </style>
        {html_table}
        """
        return styled_html

    def save_html_file(self, html_content):
        html_path = self.protocol._getExtraPath("results.html")
        with open(html_path, 'w') as f:
            f.write(html_content)

        webbrowser.open('file://' + os.path.realpath(html_path))

    
