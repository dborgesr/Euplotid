#!/bin/bash

export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket

sleep 1
cd /usr/src/app/ && node resin-wifi-connect/src/app.js --clear=false

cd /root/Euplotid
jupyter notebook --allow-root --port=80 --no-browser --ip=0.0.0.0
