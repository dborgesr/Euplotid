import os
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.figure_factory as ff
import plotly.graph_objs as go
import pandas as pd
import subprocess
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import create_engine, Date, cast
from datetime import datetime as dt
from crontab import CronTab

# Set up Dash app and database
server = Flask('applotid')
server.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////root/Euplotid/euploDB.db'
server.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(server)
#create tables, may need to make sure to not create if they exist already
db.create_all()
#start Dash app
app = dash.Dash('applotid-front', server=server,csrf_protect=False)

#Add css style to dash app
app.css.append_css({
    "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
})

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
        
def draw_ph_graph(results):
    ph_plot = go.Scatter(
        x = results["date_time"],
        y = results["ph"],
        name = 'pH',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "pH", range = [1,12])
    return fig

def draw_ec_graph(results):
    ec_plot = go.Scatter(
        x = results["date_time"],
        y = results["ec"],
        name = 'Electrical Conductivity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "MicroSiemens", range = [5,200])
    return fig

def draw_tds_plot(results):
    tds_plot = go.Scatter(
        x = results["date_time"],
        y = results["tds"],
        name = 'Total Dissolved Solids',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "Concentration", range = [0,100])

def draw_sal_plot(results):
    sal_plot = go.Scatter(
        x = results["date_time"],
        y = results["sal"],
        name = 'Salinity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "Practical salinity unit", range = [0,100])
    return fig

def draw_sg_plot(results):
    sg_plot = go.Scatter(
        x = results["date_time"],
        y = results["sg"],
        name = 'Density',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(range = [0,100]) 
    return fig

def draw_temp_plot(results):
    temp_plot = go.Scatter(
        x = results["date_time"],
        y = results["temperature"],
        name = 'Temperature',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "Farenheit", range = [-10,100])
    return fig

def draw_pressure_plot(results):
    pres_plot = go.Scatter(
        x = results["date_time"],
        y = results["pressure"],
        name = 'Pressure',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "Millibars", range = [0,1100])
    return fig

def draw_humidity_graph(results):
    hum_plot = go.Scatter(
        x = results["date_time"],
        y = results["humidity"],
        name = 'Humidity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    fig['layout']['yaxis'].update(title = "Relative Humidity", range = [0,100])
    return fig
        
def add_to_cron(minute, hour, day, month, day_week, cron_job):
    cron_euplo = CronTab(user='root')
    if cron_job:
        job  = cron_euplo.new(command=cron_job)
        job.setall(minute, hour, day, month, day_week)
        cron_euplo.write()
    df = pd.DataFrame(columns=["minute", "hour", "day", "month", "day of week", "command"])
    for saved_job in cron_euplo:
        a = saved_job.split("\t")
        job_row = [a[0],a[1],a[2],a[3],a[4],"\t".join(a[5:(len(a)-1)])]
        df.append(job_row)
    
    return df
    
def get_tentacle_readings(start_date, end_date):
    #fetch all results within a date range from DB and return dataframe
    query = tentacleRead.query.filter_by(tentacleRead.date_time >= start_date 
                                        and tentacleRead.date_time <= end_date)
    results = pd.read_sql(query.statement, query.session.bind)
    return results

def get_pisense_readings(start_date, end_date):
    #fetch all results within a date range from DB and return dataframe
    query = piSenseRead.query.filter_by(piSenseRead.date_time >= start_date 
                                        and piSenseRead.date_time <= end_date)
    results = pd.read_sql(query.statement, query.session.bind)
    return results

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
        html.H1('Euplotid Dashboard', style={'text-align': 'center'}),
        dcc.Link('Jupyter link', href='/jupyter/')
    ]),
    #Control RF connected devices
    html.Div([
        html.Div([
            html.H2('Control X10 devices')
        ]),
        html.Div([
            html.Button('RF1 ON', id='rf_on_button1'),
            html.Button('RF1 OFF', id='rf_off_button1'),
            html.P(id='placeholder1'),
            html.P(id='placeholder2')
        ], className='row'),
        html.Div([
            html.Button('RF2 ON', id='rf_on_button2'),
            html.Button('RF2 OFF', id='rf_off_button2'),
            html.P(id='placeholder3'),
            html.P(id='placeholder4')
        ], className='row'),
        html.Div([
            html.Button('RF3 ON', id='rf_on_button3'),
            html.Button('RF3 OFF', id='rf_off_button3'),
            html.P(id='placeholder5'),
            html.P(id='placeholder6')
        ], className='row'),
        html.Div([
            html.Button('RF4 ON', id='rf_on_button4'),
            html.Button('RF4 OFF', id='rf_off_button4'),
            html.P(id='placeholder7'),
            html.P(id='placeholder8')
        ], className='row'),
        html.Div([
            html.Button('RF5 ON', id='rf_on_button5'),
            html.Button('RF5 OFF', id='rf_off_button5'),
            html.P(id='placeholder9'),
            html.P(id='placeholder10')
        ], className='row')
    ]),
    
    
    #Environmental Monitoring
    html.Div([
        
        html.Div([
            html.H2('Environmental Monitoring', style={'text-align': 'center'})
        ]),
        
        html.Div([
            html.H3('Pick dates to query environment'),
            #pick dates for environmental monitoring
            dcc.DatePickerRange(
                id='date-picker-range',
                start_date=dt(2018, 1, 1),
                end_date_placeholder_text='Select a date!',
                end_date=dt(2019,1,1),
                min_date_allowed=dt(2018,1,1),
                max_date_allowed=dt(2020,1,1)
            ),
            html.Button(id='date-enviro-submit', n_clicks=0, children='Submit')
        ], className="row"),
        
        html.Div([
            html.Div([
                html.H3('pH Graph'),
                dcc.Graph(id='ph-graph')
            ], className="six columns"),
            html.Div([
                html.H3('Electrical Conductivity Graph'),
                dcc.Graph(id='ec-graph')
            ], className="six columns"),
        ], className="row"),
        
        html.Div([
            html.Div([
                html.H3('Total Dissolved Solids Graph'),
                dcc.Graph(id='tds-graph')
            ], className="six columns"),
            html.Div([
                html.H3('Salinity Graph'),
                dcc.Graph(id='sal-graph')
            ], className="six columns"),
        ], className="row"),
        
        html.Div([
            html.Div([
                html.H3('Specific Gravity Graph'),
                dcc.Graph(id='sg-graph')
            ], className="six columns"),
            html.Div([
                html.H3('Temperature Graph'),
                dcc.Graph(id='temp-graph')
            ], className="six columns"),
        ], className="row"),
        
        html.Div([
            html.Div([
                html.H3('Pressure Graph'),
                dcc.Graph(id='pressure-graph')
            ], className="six columns"),
            html.Div([
                html.H3('Humidity Graph'),
                dcc.Graph(id='humidity-graph')
            ], className="six columns")
        ], className="row")
    ]),
    
    #Cron job input and scheduling
    html.Div([
        html.H2('Schedule jobs using CRON', style={'text-align': 'center'}),
        html.Div([
            dcc.Dropdown(
                id="min-slider",
                options=[{ 'label': minute, 'value': minute } for minute in range(60)],
                value='*'
            ),
            dcc.Dropdown(
                id="hour-slider",
                options=[{ 'label': hour, 'value': hour } for hour in range(24)],
                value='*'
            ),
            dcc.Dropdown(
                id="day-slider",
                options=[{ 'label': day, 'value': day } for day in range(1,32)],
                value='*'
            ),
            dcc.Dropdown(
                id="month-picker",
                options=[
                    {'label': 'January', 'value': 'jan'},
                    {'label': 'February', 'value': 'feb'},
                    {'label': 'March', 'value': 'SF'},
                    {'label': 'April', 'value': 'apr'},
                    {'label': 'May', 'value': 'may'},
                    {'label': 'June', 'value': 'jun'},
                    {'label': 'July', 'value': 'jul'},
                    {'label': 'August', 'value': 'aug'},
                    {'label': 'September', 'value': 'sept'},
                    {'label': 'October', 'value': 'oct'},
                    {'label': 'November', 'value': 'nov'},
                    {'label': 'December', 'value': 'dec'}
                ],
                value='*'
            ),
            dcc.Dropdown(
                id="day-picker",
                options=[
                    {'label': 'Monday', 'value': 'mon'},
                    {'label': 'Tuesday', 'value': 'tue'},
                    {'label': 'Wednesday', 'value': 'wed'},
                    {'label': 'Thursday', 'value': 'thr'},
                    {'label': 'Friday', 'value': 'fri'},
                    {'label': 'Saturday', 'value': 'sat'},
                    {'label': 'Sunday', 'value': 'sun'}
                ],
                value='*'
            ),
            dcc.Input(
                id='cron-input',
                placeholder='Enter a valid command',
                type='text',
                value=''),
            html.Button(id='cron-submit', n_clicks=0, children='Submit')
            
        ],className="seven columns"),
        
        html.Div([
            html.Table(id='cron-jobs')
        ], className="six columns")
    ]),
    
])


#############################################
# Interaction Between Components / Controller
#############################################

# Update cron jobs Table
@app.callback(
    Output(component_id='cron-jobs', component_property='children'),
    [
        Input(component_id='cron-submit', component_property='n_clicks')
    ],
    [
        State('min-slider', 'value'),
        State('hour-slider', 'value'),
        State('day-slider', 'value'),
        State('day-picker', 'value'),
        State('month-picker', 'value'),
        State('cron-input', 'value')
    ]
)
def update_cron_jobs(cron_submit, minute, hour, day, month, day_week, cron_job):
    results = add_to_cron(minute, hour, day, month, day_week, cron_job)
    return generate_table(results, max_rows=50)

# Update environmental graphs
@app.callback(
    Output(component_id='ph-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_ph_graph(date_enviro_submit, start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_ph_graph(results)
    return figure

@app.callback(
    Output(component_id='ec-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_ec_graph(date_enviro_submit, start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_ec_graph(results)
    return figure

@app.callback(
    Output(component_id='tds-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_tds_graph(date_enviro_submit, start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_tds_graph(results)
    return figure

@app.callback(
    Output(component_id='sal-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_sal_graph(date_enviro_submit, start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_sal_graph(results)
    return figure

@app.callback(
    Output(component_id='sg-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_sg_graph(date_enviro_submit, start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_sg_graph(results)
    return figure

@app.callback(
    Output(component_id='temp-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_temp_graph(date_enviro_submit, start_date, end_date):
    results = get_pisense_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_temp_graph(results)
    return figure

@app.callback(
    Output(component_id='pressure-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_pressure_graph(date_enviro_submit, start_date, end_date):
    results = get_pisense_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_pressure_graph(results)
    return figure

@app.callback(
    Output(component_id='humidity-graph', component_property='figure'),
    [
        Input(component_id='date-enviro-submit', component_property='n_clicks')

    ],
    [
        State(component_id='date-picker-range', component_property='start_date'),
        State(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_humidity_graph(date_enviro_submit, start_date, end_date):
    results = get_pisense_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_humidity_graph(results)
    return figure

#Buttons for RF control
#ON
@app.callback(
    Output(component_id='placeholder1', component_property='children'),
    [Input(component_id='rf_on_button1', component_property='n_clicks')])
def rf1_on(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5264691 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder2', component_property='children'),
    [Input(component_id='rf_on_button2', component_property='n_clicks')])
def rf2_on(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5264835 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder3', component_property='children'),
    [Input(component_id='rf_on_button3', component_property='n_clicks')])
def rf3_on(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5265155 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder4', component_property='children'),
    [Input(component_id='rf_on_button4', component_property='n_clicks')])
def rf4_on(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5266691 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder5', component_property='children'),
    [Input(component_id='rf_on_button5', component_property='n_clicks')])
def rf5_on(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5272835 -l 192 -p 0", shell=True)  
    return
#OFF
@app.callback(
    Output(component_id='placeholder6', component_property='children'),
    [Input(component_id='rf_off_button1', component_property='n_clicks')])
def rf1_off(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5264700 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder7', component_property='children'),
    [Input(component_id='rf_off_button2', component_property='n_clicks')])
def rf2_off(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5264844 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder8', component_property='children'),
    [Input(component_id='rf_off_button3', component_property='n_clicks')])
def rf3_off(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5265164 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder9', component_property='children'),
    [Input(component_id='rf_off_button4', component_property='n_clicks')])
def rf4_off(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5266700 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder10', component_property='children'),
    [Input(component_id='rf_off_button5', component_property='n_clicks')])
def rf5_off(n_clicks):
    subprocess.check_output("/var/www/rfoutlet/codesend 5272844 -l 192 -p 0", shell=True)  
    return

if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0',
        port=8050
    )
