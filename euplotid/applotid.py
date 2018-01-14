import os
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.figure_factory as ff
import plotly.graph_objs as go
import pandas as pd
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import create_engine

###########################
# Data Manipulation / Model
###########################

class piSenseRead(db.Model):
    id = db.Column(db.Integer, primary_key = True)
    date_time = db.Column(db.DateTime)
    pressure = db.Column(db.Float)
    temperature = db.Column(db.Float)
    humidity = db.Column(db.Float)
    
    def __init__(self, id, date_time, pressure, temperature, humidity):
        self.id = id
        self.date_time = date_time
        self.pressure = pressure
        self.temperature = temperature
        self.humidity = humidity

class tentacleRead(db.Model):
    id = db.Column(db.Integer, primary_key = True)
    date_time = db.Column(db.DateTime)
    pH = db.Column(db.Float)
    ec = db.Column(db.Float)
    tds = db.Column(db.Float)
    sal = db.Column(db.Float)
    sg = db.Column(db.Float)
    
    def __init__(self, id, date_time, pH, ec, tds, sal, sg):
        self.id = id
        self.date_time = date_time
        self.pH = pH
        self.ec = ec
        self.tds = tds
        self.sal = sal
        self.sg = sg
        
def draw_tentacle_graph(results):
    dates = results['date']
    points = results['points'].cumsum()

    figure = go.Figure(
        data=[
            go.Scatter(x=dates, y=points, mode='lines+markers')
        ],
        layout=go.Layout(
            title='Points Accumulation',
            showlegend=False
        )
    )

    return figure
        
    
def draw_pisense_graph(results):
    dates = results['date']
    points = results['points'].cumsum()

    figure = go.Figure(
        data=[
            go.Scatter(x=dates, y=points, mode='lines+markers')
        ],
        layout=go.Layout(
            title='Points Accumulation',
            showlegend=False
        )
    )

    return figure
        
def update_cron_jobs(cron_input):
    return results
    
def get_tentacle_readings(date_input):
    return results

def get_pisense_readings(date_input):
    return results

# Set up Dash app and database
app = dash.Dash()
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////root/Euplotid/euploDB.db'
db = SQLAlchemy(app)

engine = create_engine("sqlite:////root/Euplotid/euploDB.db")  # Access the DB Engine
if not engine.dialect.has_table(engine, Variable_tableName):  # If table don't exist, Create.
    db.create_all()

app.css.append_css({
    "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
})

#########################
# Dashboard Layout / View
#########################

def generate_table(dataframe, max_rows=10):
    '''Given dataframe, return template generated using Dash components
    '''
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))]
    )

app.layout = html.Div([
    
    html.Div([
        html.H1('Euplotid Dashboard')
    ]),
    
    html.Div([
        #date selector
        dcc.Dropdown(
            id='date-picker',
            options=unique_dates,
            multi=True,
            value=[]
        ),
        # pisense output Graph
        dcc.Graph(id='pisense-graph'),
        # graph
        dcc.Graph(id='tentacle-graph')
    ], className='three columns')
    
    html.Div([
        dcc.Input(
            id='cron-input',
            placeholder='Enter a valid cronjob',
            type='text',
            value=''),
        html.Div([
            html.Table(id='cron-jobs')
        ], className='six columns')
    ],className='two columns')
])


#############################################
# Interaction Between Components / Controller
#############################################

# Update cron jobs Table
@app.callback(
    Output(component_id='cron-jobs', component_property='children'),
    [
        Input(component_id='cron-input', component_property='value')
    ]
)
def update_cron_jobs(cron_input):
    results = add_to_cron(cron_input)
    return generate_table(results, max_rows=50)


# Update piSense graph
@app.callback(
    Output(component_id='pisense-graph', component_property='figure'),
    [
        Input(component_id='date-picker', component_property='value')
    ]
)
def load_pisense_graph(dates_picked):
    results = get_pisense_readings(dates_picked)

    figure = []
    if len(results) > 0:
        figure = draw_pisense_graph(results)

    return figure

# Update tentacle graph
@app.callback(
    Output(component_id='tentacle-graph', component_property='figure'),
    [
        Input(component_id='date-picker', component_property='value')
    ]
)
def load_tentacle_graph(dates_picked):
    results = get_tentacle_readings(dates_picked)

    figure = []
    if len(results) > 0:
        figure = draw_tentacle_graph(results)

    return figure



if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0',
        port=8050
    )
