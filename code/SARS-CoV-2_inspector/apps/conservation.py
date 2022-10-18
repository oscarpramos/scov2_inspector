### MODULES IMPORT ###
import dash 
from dash.dependencies import Input, Output
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
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
                
                
# Plot results
def PlotLines(dff):
	fig = px.scatter(dff, x="PosNbr", y="Conservation", color="Annotation(s)", labels={'PosNbr':'Position in the genome', 'Conservation':'Conservation frequency (0-1): ' }, hover_data=['PosNbr','Conservation', 'Annotation(s)', 'NTmutation(s)', 'AAmutation(s)' ], title="Conservation by position (Database date: %s) <br>Total number of genomes = %i" % (InfoDict['GDD'][1], int(InfoDict['NFG'][1])))
	fig.update_layout(hoverlabel_align = 'left', height = 700)
	fig.update_yaxes(range=[-0.1, 1.1])

	return fig

### Initialization / RUN ###
# Read conservation table and comments
try:
        ConservTableFN = 'ConservationTable.tsv'
        InfoDict = getComments(ConservTableFN)
        skiplines = len(InfoDict.keys()) + 2
        extension = ConservTableFN.split(".")[-1]
        if  extension == "gz" or extension == "gzip":
                ConservDF = pd.read_csv(ConservTableFN, sep='\t', compression='gzip', skiprows=skiplines)
        else:
                ConservDF = pd.read_csv(ConservTableFN, sep='\t', skiprows=skiplines)
except FileNotFoundError:
        print("The file %s was not found !" % (filename))


NbrOfListGen = ConservDF.shape[0]

### dash output ###
layout = html.Div([
    html.H1(children='SARS-CoV-2 genome data inspector'),
    html.Br(),
    html.H2("{} : {}".format(InfoDict['GDD'][0], InfoDict['GDD'][1]), style = {'text-align': 'left'}),
    html.H5("{} : {:,}".format(InfoDict['NoG'][0], int(InfoDict['NoG'][1])).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['AUF'][0], InfoDict['AUF'][1].replace(",", ", ")), style = {'text-align': 'left'}),
    html.H5("{} : {:,}".format(InfoDict['NFG'][0], int(InfoDict['NFG'][1])).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{:,} genomic positions are listed".format(NbrOfListGen).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['RID'][0], InfoDict['RID'][1]), style = {'text-align': 'left'}),
    html.Br(),
    dash_table.DataTable(
        id='datatable-interactivity_4',
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in ConservDF.columns
        ],
         style_data={
        'whiteSpace': 'normal',
    	},
#    	css=[{
#        'selector': '.dash-spreadsheet td div',
#        'rule': '''
##            line-height: 15px;
##            max-height: 30px; min-height: 30px; height: 30px;
#	    min-width: 180px;
#            display: block;
#            overflow-y: hidden;
##        '''
#        }],
        style_header={
        'backgroundColor': 'rgb(230, 230, 230)',
        'fontWeight': 'bold'
        },
        style_cell={
        'height': 'auto',
        # all three widths are needed
        'maxWidth': '300px',
        'minWidth': '50px',
        'whiteSpace': 'normal',
	 'textOverflow': 'ellipsis',
        },
        data=ConservDF.to_dict('records'),
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
#    	 export_headers='display',
    	merge_duplicate_headers=True,
    ),
    html.Div(id='datatable-interactivity-container4B'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container4C'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container4D'),
])


@app.callback(
    Output('datatable-interactivity_4', 'style_data_conditional2'),
    Input('datatable-interactivity_4', 'selected_columns')
)

def update_styles(selected_columns):
    return [{
        'if': { 'column_id': i },
        'background_color': '#D2F3FF'
    } for i in selected_columns]

@app.callback(
	Output('datatable-interactivity-container4D', "children"),
	Input('datatable-interactivity_4', "derived_virtual_data"),
	Input('datatable-interactivity_4', "derived_viewport_selected_rows"),
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
    
    dff5 = ConservDF if rows is None else pd.DataFrame(rows)
    
    Fig = PlotLines(dff5)

    colors = ['red' if i in derived_virtual_selected_rows else '#0074D9'
              for i in range(len(dff5))]
    html.Br()
    return [
        dcc.Graph(
            id="Conservation",
            figure=Fig
        )
    ]
