import dash_html_components as html
import dash_bootstrap_components as dbc

# change to app.layout if running as single page app instead
layout = html.Div([
    dbc.Container([
        dbc.Row
            ([
            	dbc.Col(html.B(html.H1("SARS-CoV-2 genome data inspector", style={'text-align':'center'})))
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children='Please choose one of the options available from the “explore” menu above '), className="mb-4")
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children='Options:'), className="mb-6")
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children=' * Individual genomes: Explore and visualize individual genomes'), className="mb-6")
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children='* SNPs: Explore and visualize SNPs and their frequency evolution'), className="mb-6")
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children=' * Genotypes: Explore and visualize non-redundant genotypes and their frequency evolution'), className="mb-6")
            ]),
        dbc.Row
            ([
            	dbc.Col(html.H4(children=' * Conservation: Visualize conservation compared to reference genome over all genomic positions'), className="mb-6")
            ]),
    ])
])
