#!/bin/bash

diff /root/Euplotid/arduino_code/euplouino.ino /data/euplouino.ino || PROGRAMMER=1
if [ "${PROGRAMMER:-}" == "1" ]; then
  cd /root/Euplotid/arduino_code/
  make upload && cp euplouino.ino /data/
  unset PROGRAMMER
fi

export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket

sleep 1
cd /usr/src/app/ && ./wifi-connect --clear=false

#Start app for Euplotid
cd /root/Euplotid/euplotid && python applotid.py &
#Start jupyter
cd /root/Euplotid
jupyter notebook --allow-root --port=8888 --no-browser --ip=0.0.0.0
