from sense_hat import SenseHat
import subprocess, re, sys, os, pandas as pd, datetime
from applotid import db

def main():
    # Query pisense
    sense = SenseHat()
    # Temperature
    t = sense.get_temperature()
    cpu_temp = str((str(subprocess.check_output("vcgencmd measure_temp", shell=True)).replace("temp=","")).replace("'C",""))
    cpu_temp = float(re.findall("\d+\.\d+", cpu_temp)[0])
    t_c = t - ((cpu_temp - t)/1.166) # calibrate according to CPU temp
    t_c = round(((t_c/5.0)*9)+32,1) #farenheit
    # Compass
    compass = sense.get_compass()
    # Humidity
    h = sense.get_humidity()
    h = round(h, 1)
    # Barometric Pressure
    p = round(sense.get_pressure(),2)
    msg = "Euplotid says: Temp = {0}'F, Press = {1}Pa, Humidity = {2}%".format(t_c,p,h)
    sense.show_message(msg, scroll_speed=0.02)
    
    #create and save entry in DB
    pisense_dbrow = piSenseRead(datetime.datetime.now(), p, t_c, h)
    db.session.add(pisense_dbrow)
    db.session.commit()
    sense.clear()
    

if __name__ == "__main__":
    main()