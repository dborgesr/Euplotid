#!/bin/bash

modprobe i2c-dev #sensehat
cd /usr/src/app/RPi_Cam_Web_Interface && service nginx restart #microscope server
#modprobe bcm2835-v4l2 #raspicam

export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket

sleep 1
cd /usr/src/app/ && node resin-wifi-connect/src/app.js --clear=true

cd /root/Euplotid
jupyter notebook --allow-root --port=80 --no-browser --ip=0.0.0.0

while true; do
	sleep 100
done