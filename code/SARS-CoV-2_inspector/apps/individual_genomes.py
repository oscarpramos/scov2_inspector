### MODULES IMPORT ###
import dash 
from dash.dependencies import Input, Output
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import pandas as pd
import gzip

from app import app

#### FUNCTIONS ###
# get comments from a file	
def getComments(FN):
	try:
		start = False
		my_splitlines = []
		InfoDict = {}
		extension = FN.split(".")[-1]
		if extension == "gz" or extension == "gzip":
			with gzip.open(FN, mode='rb') as f:
				for line in f:
					line = line.decode("utf-8")
					if line.split("|")[0][1:4] == 'END':	
						start = False
						break
					if start == True:
						my_splitlines.append(line)	
					if line.split("|")[0][1:6] == 'START':
						start = True
					
			f.close()
		else:
			with open(FN, mode='r') as f:
				for line in f:
					if line.split("|")[0][1:4] == 'END':	
						start = False
						break
					if start == True:
						my_splitlines.append(line)	
					if line.split("|")[0][1:6] == 'START':
						start = True
			f.close()

		for comment in my_splitlines:
			code = comment.split("|")[0][1:4]
			info = comment.split(":")[-1].strip()
			desc = comment.split("|")[-1].split(":")[0]
			InfoDict[code] = [desc, info]

		return InfoDict
		
	except FileNotFoundError:
                print("The file %s was not found !" % (FN))

# Reformat a mutation list to make it fit in the hover panel
def reformatMutaList(genotype, NbrOfMutaLine):
    MutaTextArray = []
    for muta in genotype:
        mutas = muta.split("|")
        a = 1
        MutaText = ""
        
        for mut in mutas:
                Concat = [MutaText, mut]
                if MutaText == "":
                        MutaText = mut
                elif a > NbrOfMutaLine:
                        MutaText = "<br>".join(Concat)
                        a = 1
                else:
                        MutaText = "|".join(Concat)
                        a += 1
        MutaTextArray.append(MutaText)
       
    MutaTextArray
    return MutaTextArray

## Initialization
# Read genomes table and comments
try:
        IndGenTableFN = 'IndividualGenomeList_100K.tsv.gz'
        InfoDict = getComments(IndGenTableFN)
        skiplines = len(InfoDict.keys()) + 2
        extension = IndGenTableFN.split(".")[-1]
        if  extension == "gz" or extension == "gzip":
                GenomeDF = pd.read_csv(IndGenTableFN, sep='\t', compression='gzip', skiprows=skiplines, low_memory=False)
        else:
                GenomeDF = pd.read_csv(IndGenTableFN, sep='\t', skiprows=skiplines, low_memory=False)
except FileNotFoundError:
        print("The file %s was not found !" % (filename))

NbrOfListGen = GenomeDF.shape[0]

### dash output ###
layout = html.Div([
    html.H1(children='SARS-CoV-2 genome data inspector'),
    html.Br(),
    html.H2("{} : {}".format(InfoDict['GDD'][0], InfoDict['GDD'][1]), style = {'text-align': 'left'}),
    html.H5("{} : {:,}".format(InfoDict['NoG'][0], int(InfoDict['NoG'][1])).replace(",", " "), style = {'text-align': 'left'}),
#    html.H5("{} : {}".format(InfoDict['AUF'][0], InfoDict['AUF'][1].replace(",", ", ")), style = {'text-align': 'left'}),              # can be uncommented if the indiv. gen. list was filtered
#    html.H5("{} : {:,}".format(InfoDict['NFG'][0], int(InfoDict['NFG'][1])).replace(",", " "), style = {'text-align': 'left'}),        # can be uncommented if the indiv. gen. list was filtered
    html.H5("{:,} genomes are listed".format(NbrOfListGen).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['RID'][0], InfoDict['RID'][1]), style = {'text-align': 'left'}),
    
    html.Br(),
    dash_table.DataTable(

        id='datatable-interactivity_1',
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in GenomeDF.columns
        ],
         style_data={
        'whiteSpace': 'normal',
    	},
#    	css=[{
#        'selector': '.dash-spreadsheet td div',
#        'rule': '''
#            line-height: 15px;
#            max-height: 400px; min-height: 60px; height: 90px;
#            display: block;
#            overflow-y: hidden;
#        '''
#        }],
        style_header={
        'backgroundColor': 'rgb(230, 230, 230)',
        'fontWeight': 'bold'
        },
        style_cell={
        'height': 'auto',
        # all three widths are needed
        'maxWidth': '300px',
        'minWidth': '100px',
        'whiteSpace': 'normal',
	 'textOverflow': 'ellipsis',
        },
        data=GenomeDF.to_dict('records'),
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        row_selectable="multi",
        row_deletable=True,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 25,
#        export_format='csv',
#     	 export_headers='display',
    	merge_duplicate_headers=True,
    ),
    html.Div(id='datatable-interactivity-container1B'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container1C'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container1D'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container1E'),
    html.Br(),
    html.Br(),
    html.Br(),
])


@app.callback(
    Output('datatable-interactivity_1', 'style_data_conditional'),
    Input('datatable-interactivity_1', 'selected_columns')
)

def update_styles(selected_columns):
    return [{
        'if': { 'column_id': i },
        'background_color': '#D2F3FF'
    } for i in selected_columns]

@app.callback(
    Output('datatable-interactivity-container1C', "children"),
    Input('datatable-interactivity_1', "derived_viewport_data"),
    Input('datatable-interactivity_1', "derived_viewport_selected_rows"),
    )
    
def update_graphs(rows, derived_virtual_selected_rows):
    # When the table is first rendered, `derived_virtual_data` and
    # `derived_virtual_selected_rows` will be `None`. This is due to an
    # idiosyncrasy in Dash (unsupplied properties are always None and Dash
    # calls the dependent callbacks when the component is first rendered).
    # So, if `rows` is `None`, then the component was just rendered
    # and its value will be the same as the component's dataframe.
    # Instead of setting `None` in here, you could also set
    # `derived_virtual_data=df.to_rows('dict')` when you initialize
    # the component.
    if derived_virtual_selected_rows is None:
        derived_virtual_selected_rows = []

    dff = GenomeDF if rows is None else pd.DataFrame(rows)

    dff['PointMutations'] = reformatMutaList(dff['PointMutations'], 11) # genotype and NbrOfMuta per line

    colors = ['red' if i in derived_virtual_selected_rows else '#0074D9'
              for i in range(len(dff))]
    html.Br()
    return [
        dcc.Graph(
            id=column,
            figure={
                "data": [
                    {
                        "x": dff["GISAID_ID"],
                        "y": dff[column],
                        "type": "scatter",
                        "mode": "markers",
                        "marker": {"color": colors},
                        "name": "Genome",
                        "hovermode": "closest",
                        "hovertemplate": "<i><b>GISAID_ID</b></i>: %{x}<br><i>Y axis:</i> %{y}<br><i>Collection date</i>: %{hovertext}<br><br><i><b>Genotype(NT mutation;AA mutation;Annotation)</b></i><br>%{customdata}<br>",
                        "hovertext": dff["Collection_Date"],
                        "customdata": dff['PointMutations'],
                    }
                ],
                "layout": {
                    "xaxis": {"automargin": True},
                    "yaxis": {
                        "automargin": True,
                        "title": {"text": column}
                    },
                    "height": 250,
                    "margin": {"t": 10, "l": 10, "r": 10},
                },
            },
        )
        # check if column exists - user may have deleted it
        # If `column.deletable=False`, then you don't
        # need to do this check.
        for column in ["#Mutations", "Collection_Date", "Country", "Clade"] if column in dff
    ]
    
