import os
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.figure_factory as ff
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import subprocess
#from picamera.array import PiRGBArray
#from picamera import PiCamera
import time
#import cv2
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import create_engine, Date, cast
from datetime import datetime as dt
import time
from crontab import CronTab

# Set up Dash app and database
server = Flask('applotid')
server.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////root/Euplotid/euploDB.db'
server.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(server)


#start webcam
#camera = PiCamera()
#camera.start_preview()
#camera.capture('t.jpg')
 
# initialize the camera and grab a reference to the raw camera capture
#camera = PiCamera()
#rawCapture = PiRGBArray(camera)

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
    data = go.Scatter(
        x = results["date_time"],
        y = results["ph"],
        name = 'pH',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='pH',
            range=[1,12]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_ec_graph(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["ec"],
        name = 'Electrical Conductivity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='MicroSiemens',
            range=[5,200]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_tds_plot(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["tds"],
        name = 'Total Dissolved Solids',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='Concentration',
            range=[0,100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig
    
def draw_sal_plot(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["sal"],
        name = 'Salinity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='Practical salinity unit',
            range=[0,100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_sg_plot(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["sg"],
        name = 'Density',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            range=[0,100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_temp_plot(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["temperature"],
        name = 'Temperature',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='Farenheit',
            range=[-10,100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_pressure_plot(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["pressure"],
        name = 'Pressure',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='Millibars',
            range=[0,1100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig

def draw_humidity_graph(results):
    data = go.Scatter(
        x = results["date_time"],
        y = results["humidity"],
        name = 'Humidity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False),
    layout = go.Layout(
        xaxis=dict(
            title='Date'
        ),
        yaxis=dict(
            title='Relative Humidity',
            range=[0,100]
        )
    )
    
    fig = go.Figure(data=data, layout=layout)
    
    return fig
        
def add_to_cron(minute, hour, day, month, day_week, cron_job):
    cron_euplo = CronTab(user='root')
    if cron_job:
        job  = cron_euplo.new(command=cron_job)
        job.setall(minute, hour, day, month, day_week)
        if job.is_valid():
            cron_euplo.write()
        else:
            job.clear()
    return
    
def clear_cron():
    cron_euplo = CronTab(user='root')
    cron_euplo.remove_all()
    cron_euplo.write()
    return

def get_cron_table():
    cron_euplo = CronTab(user='root')
    df = pd.DataFrame(columns=["minute", "hour", "day", "month", "day of week", "command"])
    job_arr = np.empty((0,6), str)
    for saved_job in cron_euplo:
        a = str(saved_job).split(" ")
        job_row = np.reshape(np.array([a[0],a[1],a[2],a[3],a[4]," ".join(a[5:(len(a))])]),(1,6))
        job_arr = np.vstack((job_arr,job_row))
    job_df = pd.DataFrame(job_arr, columns=["minute", "hour", "day", "month", "day of week", "command"])
    
    return job_df

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

#Cron job input and scheduling options
min_opts = [dict(label=str(minute), value=str(minute)) for minute in range(60)]
min_opts.append(dict(label=str('*'), value=str('*')))
hr_opts = [dict(label=str(hour), value=str(hour)) for hour in range(24)]
hr_opts.append(dict(label=str('*'), value=str('*')))
dat_opts = [dict(label=str(day), value=str(day)) for day in range(1,32)]
dat_opts.append(dict(label=str('*'), value=str('*')))

rf_code_button = ['RF1','RF2','RF3','RF4','RF5','NA','NA','NA','NA','NA']
rf_code_on = [5264691,5264835,5265155,5266691,5272835,267571,267715,268035,269571,275715]
rf_code_on_time = [192,192,192,192,192,185,187,186,186,186]
rf_code_off = [5264700,5264844,5265164,5266700,5272844,267580,267724,268044,269580,275724]
rf_code_off_time = [192,192,192,192,192,185,187,186,186,186]
columns_rf = { "RF Button":rf_code_button, "RF code ON":rf_code_on, "RF code ON timing":rf_code_on_time, "RF code OFF":rf_code_off, "RF code OFF timing":rf_code_off_time }
rf_df = pd.DataFrame(data=columns_rf)

app.layout = html.Div([
    
    html.Div(
            [
                html.H1(
                    'Euplotid Dashboard',
                    className='two columns',
                ),
                html.Img(
                    src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/Title_slide.png",
                    className='one columns',
                    style={
                        'height': '200',
                        'width': '200',
                        'float': 'right',
                        'position': 'relative',
                    },
                ),
                dcc.Link('Jupyter link', href='/jupyter/'),
                html.Div([
                    html.H3('Current time', id='date_now'),
                    dcc.Interval(
                        id='interval-component',
                        interval=1*1000, # in milliseconds
                        n_intervals=0
                    )
                ], className='seven columns'),
            ],
            className='row'),
    
        
    #Control RF connected devices
    html.Div([
        
        html.H2('Control X10 devices', className='row'),
        #RF codes for known X10 devices
        html.Div([
            generate_table(rf_df,max_rows=50)
        ], className="five columns"),
        
        #Buttons for RF devices
        html.Div([
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
        ], className='two columns'),
    ],className='row'),
    
    html.Div([
        html.H2('Schedule jobs using CRON', style={'text-align': 'center'}, className='row'),
        html.Div([
            dcc.Dropdown(
                id="min-slider",
                options=min_opts,
                value='*'
            ),
            dcc.Dropdown(
                id="hour-slider",
                options=hr_opts,
                value='*'
            ),
            dcc.Dropdown(
                id="day-slider",
                options=dat_opts,
                value='*'
            ),
            dcc.Dropdown(
                id="month-picker",
                options=[
                    {'label': 'January', 'value': '1'},
                    {'label': 'February', 'value': '2'},
                    {'label': 'March', 'value': '3'},
                    {'label': 'April', 'value': '4'},
                    {'label': 'May', 'value': '5'},
                    {'label': 'June', 'value': '6'},
                    {'label': 'July', 'value': '7'},
                    {'label': 'August', 'value': '8'},
                    {'label': 'September', 'value': '9'},
                    {'label': 'October', 'value': '10'},
                    {'label': 'November', 'value': '11'},
                    {'label': 'December', 'value': '12'},
                    {'label': '*', 'value': '*'}
                ],
                value='*'
            ),
            dcc.Dropdown(
                id="day-picker",
                options=[
                    {'label': 'Monday', 'value': '0'},
                    {'label': 'Tuesday', 'value': '1'},
                    {'label': 'Wednesday', 'value': '2'},
                    {'label': 'Thursday', 'value': '3'},
                    {'label': 'Friday', 'value': '4'},
                    {'label': 'Saturday', 'value': '5'},
                    {'label': 'Sunday', 'value': '6'},
                    {'label': '*', 'value': '*'}
                ],
                value='*'
            ),
            dcc.Input(
                id='cron-input',
                placeholder='Enter a valid command',
                type='text',
                value=''),
            html.Button(id='cron-submit', n_clicks=0, children='Submit'),
            html.P(id='placeholdercron2')
            
        ],className="seven columns"),
        
        html.Div([
            html.Table(id='cron-jobs')
        ], className="six columns"),
        html.Div([
            html.Button('Erase all CRON jobs', id='cron-erase'),
            html.P(id='placeholdercron'),
            html.Button('Refresh CRON table', id='cron-refresh')
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
    ])
    
])

#############################################
# Interaction Between Components / Controller
#############################################

# Update cron jobs Table
@app.callback(
    Output(component_id='placeholdercron2', component_property='children'),
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
	if (cron_submit):
    	results = add_to_cron(minute, hour, day, month, day_week, cron_job)
    return 

# Erase CRON jobs
@app.callback(
    Output(component_id='placeholdercron', component_property='children'),
    [
        Input(component_id='cron-erase', component_property='n_clicks')
    ]
)
def erase_cron_jobs(n_clicks):
	if (n_clicks):
    	clear_cron()
    return
    
# Refresh CRON table
@app.callback(
    Output(component_id='cron-jobs', component_property='children'),
    [
        Input(component_id='cron-refresh', component_property='n_clicks'),
        Input(component_id='cron-submit', component_property='n_clicks'),
        Input(component_id='cron-erase', component_property='n_clicks')
    ]
)
def refresh_cron_table(cron_refresh_clicks, cron_submit_clicks, cron_erase_clicks):
	if (cron_refresh_clicks) and (cron_submit_clicks) and (cron_erase_clicks):
    	results = get_cron_table()
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
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5264691 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder2', component_property='children'),
    [Input(component_id='rf_on_button2', component_property='n_clicks')])
def rf2_on(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5264835 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder3', component_property='children'),
    [Input(component_id='rf_on_button3', component_property='n_clicks')])
def rf3_on(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5265155 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder4', component_property='children'),
    [Input(component_id='rf_on_button4', component_property='n_clicks')])
def rf4_on(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5266691 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder5', component_property='children'),
    [Input(component_id='rf_on_button5', component_property='n_clicks')])
def rf5_on(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5272835 -l 192 -p 0", shell=True)  
    return
#OFF
@app.callback(
    Output(component_id='placeholder6', component_property='children'),
    [Input(component_id='rf_off_button1', component_property='n_clicks')])
def rf1_off(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5264700 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder7', component_property='children'),
    [Input(component_id='rf_off_button2', component_property='n_clicks')])
def rf2_off(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5264844 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder8', component_property='children'),
    [Input(component_id='rf_off_button3', component_property='n_clicks')])
def rf3_off(n_clicks):
    if (n_clicks):
        (n_clicks) and subprocess.check_output("/var/www/rfoutlet/codesend 5265164 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder9', component_property='children'),
    [Input(component_id='rf_off_button4', component_property='n_clicks')])
def rf4_off(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5266700 -l 192 -p 0", shell=True)  
    return
@app.callback(
    Output(component_id='placeholder10', component_property='children'),
    [Input(component_id='rf_off_button5', component_property='n_clicks')])
def rf5_off(n_clicks):
    if (n_clicks):
        subprocess.check_output("/var/www/rfoutlet/codesend 5272844 -l 192 -p 0", shell=True)  
    return
    
@app.callback(Output('date_now', 'children'),
              [Input('interval-component', 'n_intervals')])
def update_time(n): 
    now = time.strftime("%c")
    ## date and time representation
    return("Current date & time " + time.strftime("%c"))

if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0',
        port=80
    )
