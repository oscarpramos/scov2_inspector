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
# transfer the content of a file to the RAM  
def get_file_content(filename):
        with open(filename) as f:
                data = f.read()
                my_splitlines = data.splitlines()
                f.close()
                return my_splitlines
                
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

# Get dictionaries of continent, countries and their codes
def get_LocDict(LocTable, scale):
	LocDict = {}
	if scale.lower() == "continent":
		for line in LocTable:
			fields = line.rsplit("\t")
			if fields[0] not in LocDict.keys():
				LocDict[fields[0]] =  fields[1]
	elif scale.lower() == "country":
		for line in LocTable:
			fields = line.rsplit("\t")
			if fields[2] not in LocDict.keys():
				LocDict[fields[2]] =  fields[4]
#			else:
#				print("This country, %s, was already seen !" % (fields[2]))
	else:
		print("Error: Scale not found ! ")
		
	return LocDict

## Initialization
# Get continent and country codes 
LocTable = get_file_content("country-and-continent-codes-list.tsv")
ContinentDict = get_LocDict(LocTable[1:], "Continent")
CountryDict = get_LocDict(LocTable[1:], "Country")
CountryCodeDict = {value : key for (key, value) in CountryDict.items()}

# Read genomes table and comments
try:
        IndGenTableFN = 'snpsTable.tsv.gz'
        InfoDict = getComments(IndGenTableFN)
        skiplines = len(InfoDict.keys()) + 2
        extension = IndGenTableFN.split(".")[-1]
        if  extension == "gz" or extension == "gzip":
                SNPsDF = pd.read_csv(IndGenTableFN, sep='\t', compression='gzip', skiprows=skiplines, dtype={'Clade': str})
        else:
                SNPsDF = pd.read_csv(IndGenTableFN, sep='\t', skiprows=skiplines, dtype={'Clade': str})
except FileNotFoundError:
        print("The file %s was not found !" % (filename))

NbrOfListGen = SNPsDF.shape[0]
    
# Plot results
def PlotLines(dff):
        # find index of time frame colmumns 
	ColumnNames = [name for name in dff.columns]
	MyIndexList = []
	for n in range(0, len(ColumnNames)):
		splitname = ColumnNames[n].split("_")
		if len(splitname) > 2 and splitname[2] == 'to':
			MyIndexList.append(n)
	FirstIndex = MyIndexList[0]
	LastIndex = MyIndexList[-1]
	
	# get number of TF and Table lines
	NumberOfTF = len(MyIndexList)
	NumberOfVars = dff.iloc[:,FirstIndex:LastIndex].shape[0]
	
	# create x axis values that correspond to each time frame
	xall = [n for n in range (1, NumberOfTF+1)] # create a list of values with the number of values = number of collumns to be analyzed (number of time frames)
	
	#construct the graph
	fig = go.Figure(data=go.Scatter(x=xall, y=dff.iloc[0:1,FirstIndex:LastIndex + 1 ].to_numpy()[0] , mode='lines+markers', name=str(dff['SNP'].values[0]), showlegend=True))
	for i in range (1, NumberOfVars):
		fig.add_trace(go.Scatter(x=xall, y=dff.iloc[i:i+1,FirstIndex:LastIndex + 1].to_numpy()[0], mode='lines+markers',  name=str(dff['SNP'].values[i]), showlegend=True))
	fig.update_layout(hoverlabel_align = 'right', title = "Frequency evolution of the displayed SNPs over time frames", title_x=0.5, xaxis_title="Time frame", yaxis_title="Frequency (0-1)", legend_title="Mutation")
	
	return fig

### dash output ###
layout = html.Div([
    html.H1(children='SARS-CoV-2 genome data inspector'),
    html.Br(),
    html.H2("{} : {}".format(InfoDict['GDD'][0], InfoDict['GDD'][1]), style = {'text-align': 'left'}),
    html.H5("{} : {:,}".format(InfoDict['NoG'][0], int(InfoDict['NoG'][1])).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['AUF'][0], InfoDict['AUF'][1].replace(",", ", ")), style = {'text-align': 'left'}),
    html.H5("{} : {:,}".format(InfoDict['NFG'][0], int(InfoDict['NFG'][1])).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{:,} emerging SNPs are listed".format(NbrOfListGen).replace(",", " "), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['RID'][0], InfoDict['RID'][1]), style = {'text-align': 'left'}),
    html.H5("{} : {}".format(InfoDict['TFD'][0], InfoDict['TFD'][1]), style = {'text-align': 'left'}),
    html.Br(),
    dash_table.DataTable(
        id='datatable-interactivity_3',
        columns=[
            {"name": i, "id": i, "deletable": True, "selectable": True} for i in SNPsDF.columns
        ],
         style_data={
        'whiteSpace': 'normal',
    	},
#    	css=[{
#        'selector': '.dash-spreadsheet td div',
#        'rule': '''
#            line-height: 15px;
#            max-height: 30px; min-height: 30px; height: 30px;
#            max-width: 300px; min-width: 90px; width: 180px;
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
        'minWidth': '240px',
        'whiteSpace': 'normal',
	 'textOverflow': 'ellipsis',
        },
        data=SNPsDF.to_dict('records'),
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
    html.Div(id='datatable-interactivity-container3B'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container3C'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div(id='datatable-interactivity-container3D'),
    html.Br(),
    html.Br(),
    html.Br(),
    html.H6(("*** Please select at least one SNP to show the related information on the map ***"), style = {'text-align': 'center'}),
    html.Div(id='datatable-interactivity-container3E'),

])


@app.callback(
    Output('datatable-interactivity_3', 'style_data_conditional2'),
    Input('datatable-interactivity_3', 'selected_columns')
)

def update_styles(selected_columns):
    return [{
        'if': { 'column_id': i },
        'background_color': '#D2F3FF'
    } for i in selected_columns]

@app.callback(
	Output('datatable-interactivity-container3C', "children"),
	Input('datatable-interactivity_3', "derived_viewport_data"),
	Input('datatable-interactivity_3', "derived_viewport_selected_rows"),
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
        
    dff3 = SNPsDF if rows is None else pd.DataFrame(rows)
    
#    dff3['FirstGenome'] = [ Vargenomes.rsplit(',')[0] for Vargenomes in dff3['GenomesList']]
    dff3['Position'] = [ Vargenomes.split(';')[0][1:-1] for Vargenomes in dff3['SNP']]
    
    colors = ['red' if i in derived_virtual_selected_rows else '#0074D9'
              for i in range(len(dff3))]
    
    
    html.Br()
    html.Br()
    return [
        dcc.Graph(
            id=column,
            figure={
                "data": [
                    {
                    	"type": "scatter",
                    	"mode": "markers",
                    	"name": "Slope of variants",
                    	"marker": {"color": colors},
                        "x": dff3["Position"],
                        "y": dff3[column],
                        "hovermode": "closest",
                        "hovertemplate": "<b>Mutation </b>: %{customdata}<br><br><i>Slope</i>: %{y:.6f}<br><i>Number of genomes</i>: %{hovertext}<br><br><i>",
                        "hovertext": dff3["#Genomes"], 
                        "customdata": dff3['SNP'],
#                        "type": "bar",
#                        "marker": {"color": "#0074D9"},
                    }
                ],
                "layout": {
                    "title": "<br>Slope of displayed SNPs <br>Positive values = emerging SNPs; Negative values = submerging SNPs ",
                    "xaxis": {"automargin": True, "title": {"text": "Position in genome"}},
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
        for column in ["Slope"] if column in dff3
    ]
    
@app.callback(
	Output('datatable-interactivity-container3D', "children"),
	Input('datatable-interactivity_3', "derived_viewport_data"),
	Input('datatable-interactivity_3', "derived_viewport_selected_rows"),
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
    
    dff4 = SNPsDF if rows is None else pd.DataFrame(rows)
    
#    dff4['SNP'] = [ dff4['SNP'][n] + " => " + dff4['Aa_mutation'][n].replace(",", " => ") for n in range(len(dff4['SNP']))]
#    dff4['FirstGenome'] = [ Vargenomes.rsplit(',')[0] for Vargenomes in dff4['GenomesList']]
    
        
    Fig = PlotLines(dff4)

    colors = ['red' if i in derived_virtual_selected_rows else '#0074D9'
              for i in range(len(dff4))]
    html.Br()
    return [
        dcc.Graph(
            id=column,
            figure=Fig
        )
        # check if column exists - user may have deleted it
        # If `column.deletable=False`, then you don't
        # need to do this check.
        for column in ["Slope"] if column in dff4
    ]
    
    
@app.callback(
	Output('datatable-interactivity-container3E', "children"),
	Input('datatable-interactivity_3', "derived_viewport_data"),
	Input('datatable-interactivity_3', "derived_viewport_selected_rows"),
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
    
    MyGraphDict = {}
    CountriesList = []
    CountList = []
    FreqList = []
    for index in derived_virtual_selected_rows:
    	CountriesData = rows[index]['Country']
    	if CountriesData != None:
    		CountriesDataList = CountriesData.rsplit(";")
	    	for CountryItem in CountriesDataList:
	    		CountryCode = CountryDict[CountryItem.rsplit(":")[0]]
	    		CountryCount = int(CountryItem.rsplit(":")[1])
	    		CountryFreq = float(CountryItem.rsplit(":")[2])
	    		if CountryCode not in MyGraphDict.keys():
	    			MyGraphDict[CountryCode] = [CountryFreq, CountryFreq, 1, CountryCount]
	    		else:
	    			MyGraphDict[CountryCode][3] = MyGraphDict[CountryCode][3] + CountryCount
	    			MyGraphDict[CountryCode][2] += 1
	    			MyGraphDict[CountryCode][1] = MyGraphDict[CountryCode][1] + CountryFreq
	    			MyGraphDict[CountryCode][0] = MyGraphDict[CountryCode][1] / MyGraphDict[CountryCode][2]
    	else:
    		pass
	   
    for country in  MyGraphDict.keys():
    	CountriesList.append(country)
    	FreqList.append(MyGraphDict[country][0])
    	CountList.append(MyGraphDict[country][3])

    MapDF = pd.DataFrame()
    MapDF['Countries'] =  CountriesList
    MapDF['Frequencies'] = FreqList
    MapDF['Count'] = CountList
    
   
#    fig = px.choropleth(data_frame = MapDF,
#                    locations= "Countries",
#                    color= "Frequencies",  # value in column 'Confirmed' determines color
#                    hover_name= "Countries",
#                    color_continuous_scale= 'RdYlGn',  #  color scale red, yellow green
#                    )
#    fig.update_layout(title = "Frequencies of mutations per country", title_x=0.5)

    fig = go.Figure(data=go.Choropleth(
    	locations=MapDF['Countries'],
    	z=MapDF['Frequencies'].astype(float),
    	colorscale='Reds',
    	autocolorscale=False,
    	text=MapDF['Count'], # hover text
    	marker_line_color='white', # line markers between states
    	colorbar_title="Frequencies",
    	
    ),
    layout = {
            "height": 700,  # px
        }
   )

    fig.update_layout(title = "Count and frequency of selected mutations in the world<br>- Across all time frames - ", title_x=0.5)

    html.Br()
    return [
        dcc.Graph(
            id="Frequencies",
            figure=fig,
            
        )
    ]
