import os
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.figure_factory as ff
import plotly.graph_objs as go
import pandas as pd
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import create_engine, Date, cast
from datetime import datetime as dt
from crontab import CronTab

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
    
    #Generate plots
    ph_plot = go.Scatter(
        x = results["date_time"],
        y = results["ph"],
        name = 'pH',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    ec_plot = go.Scatter(
        x = results["date_time"],
        y = results["ec"],
        name = 'Electrical Conductivity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    tds_plot = go.Scatter(
        x = results["date_time"],
        y = results["tds"],
        name = 'Total Dissolved Solids',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    sal_plot = go.Scatter(
        x = results["date_time"],
        y = results["sal"],
        name = 'Salinity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    sg_plot = go.Scatter(
        x = results["date_time"],
        y = results["sg"],
        name = 'Density',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    
    fig = tools.make_subplots(rows=2, cols=3)
                              
    fig.append_trace(ph_plot, 1, 1)
    fig.append_trace(ec_plot, 1, 2)
    fig.append_trace(sal_plot, 1, 3)
    fig.append_trace(sg_plot, 2, 1)
    fig.append_trace(tds_plot, 2, 2)
                              
    fig['layout']['yaxis1'].update(title = "pH", range = [1,12])
    fig['layout']['yaxis2'].update(title = "MicroSiemens", range = [5,200])
    fig['layout']['yaxis3'].update(title = "Concentration", range = [0,100])
    fig['layout']['yaxis4'].update(title = "Density", range = [0,100])
    fig['layout']['yaxis5'].update(title = "Parts per Million", range = [0,100]) 
    
    fig['layout'].update(height=700, width=1000, title='Tentacle Readings')
    
    return fig
        
    
def draw_pisense_graph(results):
    
    #Generate plots
    temp_plot = go.Scatter(
        x = results["date_time"],
        y = results["temperature"],
        name = 'Temperature',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    pres_plot = go.Scatter(
        x = results["date_time"],
        y = results["pressure"],
        name = 'Pressure',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
    hum_plot = go.Scatter(
        x = results["date_time"],
        y = results["humidity"],
        name = 'Humidity',
        line = dict(
            color = ('rgb(22, 96, 167)'),
            width = 4,),
        showlegend = False)
                              
    fig = tools.make_subplots(rows=1, cols=3)
                              
    fig.append_trace(temp_plot, 1, 1)
    fig.append_trace(pres_plot, 1, 2)
    fig.append_trace(hum_plot, 1, 3)
                              
    fig['layout']['yaxis1'].update(title = "Farenheit", range = [-10,100])
    fig['layout']['yaxis2'].update(title = "Millibars", range = [0,1100])
    fig['layout']['yaxis3'].update(title = "Relative Humidity", range = [0,100])
    
    fig['layout'].update(height=700, width=1000, title='PiSense Readings')
    
    return fig
        
def add_to_cron(minute, hour, day, month, day_week, cron_job):
    cron_euplo = CronTab(user="Euplotid")
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
    
    #Environmental monitoring
    html.Div([
        #date selector
        dcc.DatePickerRange(
            id='date-picker-range',
            start_date=dt(2018, 1, 1),
            end_date_placeholder_text='Select a date!',
            end_date=dt(2019,1,1),
            min_date_allowed=dt(2018,1,1),
            max_date_allowed=dt(2020,1,1)
        ),
        # pisense output Graph
        dcc.Graph(id='pisense-graph'),
        # graph
        dcc.Graph(id='tentacle-graph')
    ], className='one column'),
    
    #Cron job input and scheduling
    html.Div([
        dcc.Slider(
            id="min-slider",
            min=0,
            max=59,
            marks={i: 'Minute {}'.format(i) for i in range(60)},
            value='-1'
        ),
        dcc.Slider(
            id="hour-slider",
            min=0,
            max=23,
            marks={i: 'Hour {}'.format(i) for i in range(23)},
            value='-1'
        ),        
        dcc.Slider(
            id="day-slider",
            min=1,
            max=31,
            marks={i: 'Day of month {}'.format(i) for i in range(32)},
            value='-1'
        ),       
        dcc.Dropdown(
            id="month-picker"
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
        html.Button('Submit', id='cron-submit'),
        html.Div([
            html.Table(id='cron-jobs')
        ], className='six columns')
    ],className='two columns'),
    
    #Control RF connected devices
    html.Div([
        html.Button('RF1 ON', id='rf_on_button1'),
        html.Button('RF1 OFF', id='rf_off_button1'),
        html.Button('RF2 ON', id='rf_on_button2'),
        html.Button('RF2 OFF', id='rf_off_button2'),
        html.Button('RF3 ON', id='rf_on_button3'),
        html.Button('RF3 OFF', id='rf_off_button3'),
        html.Button('RF4 ON', id='rf_on_button4'),
        html.Button('RF4 OFF', id='rf_off_button4'),
        html.Button('RF5 ON', id='rf_on_button5'),
        html.Button('RF5 OFF', id='rf_off_button5')
    ], className='two columns')
])


#############################################
# Interaction Between Components / Controller
#############################################

# Update cron jobs Table
@app.callback(
    Output(component_id='cron-jobs', component_property='children'),
    [
        Input(component_id='cron-submit', component_property='value')
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

# Update piSense graph
@app.callback(
    Output(component_id='pisense-graph', component_property='figure'),
    [
        Input(component_id='date-picker-range', component_property='start_date'),
        Input(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_pisense_graph(start_date, end_date):
    results = get_pisense_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_pisense_graph(results)
    return figure

# Update tentacle graph
@app.callback(
    Output(component_id='tentacle-graph', component_property='figure'),
    [
        Input(component_id='date-picker-range', component_property='start_date'),
        Input(component_id='date-picker-range', component_property='end_date')
    ]
)
def load_tentacle_graph(start_date, end_date):
    results = get_tentacle_readings(start_date, end_date)
    figure = []
    if len(results) > 0:
        figure = draw_tentacle_graph(results)
    return figure

#Buttons for RF control
#ON
@app.callback( ,Input(component_id='rf_on_button1', 'n_clicks'))
def rf1_on(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5264691 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_on_button2', 'n_clicks'))
def rf2_on(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5264835 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_on_button3', 'n_clicks'))
def rf3_on(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5265155 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_on_button4', 'n_clicks'))
def rf4_on(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5266691 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_on_button5', 'n_clicks'))
def rf5_on(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5272835 -l 192 -p 0", shell=True)  
    return
#OFF
@app.callback( ,Input(component_id='rf_off_button1', 'n_clicks'))
def rf1_off(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5264700 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_off_button2', 'n_clicks'))
def rf2_off(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5264844 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_off_button3', 'n_clicks'))
def rf3_off(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5265164 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_off_button4', 'n_clicks'))
def rf4_off(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5266700 -l 192 -p 0", shell=True)  
    return
@app.callback( ,Input(component_id='rf_off_button5', 'n_clicks'))
def rf5_off(n_clicks):
    subprocess.check_output("var/www/rfoutlet/codesend 5272844 -l 192 -p 0", shell=True)  
    return

if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0',
        port=8050
    )
