#!/bin/bash

#diff /root/Euplotid/arduino_code/euplouino.ino /data/euplouino.ino || PROGRAMMER=1
if [ "${EUPLO_ARDUINO:-}" == "true" ] ; then
  cd /root/Euplotid/arduino_code/
  make upload && cp euplouino.ino /data/
  unset PROGRAMMER
fi

if [ "${EUPLO_RPICAM:-}" == "true" ]; then
	cd /usr/src/app/RPi_Cam_Web_Interface && ./start.sh
fi


export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket

sleep 1
#cd /usr/src/app/ && ./wifi-connect --clear=false

if [ "${EUPLO_DEPLOY:-}" == "true" ]; 
	then /venv/bin/python3 /app/applotid.py &
fi

if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" != "nanotid" ]; 
	then jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0
	else /bin/ash
fi