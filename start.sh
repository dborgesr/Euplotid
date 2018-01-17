#!/bin/bash

diff /Euplotid/arduino_code/euplouino.ino /data/euplouino.ino || PROGRAMMER=1
if [ "${PROGRAMMER:-}" == "1" ]; then
  pushd /Euplotid/arduino_code/
  make upload && cp euplouino.ino /data/
  unset PROGRAMMER
  popd
fi

export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket

sleep 1
cd /usr/src/app/ && node resin-wifi-connect/src/app.js --clear=false
cd /root/Euplotid/euplotid && python applotid.py &

cd /root/Euplotid
jupyter notebook --allow-root --port=80 --no-browser --ip=0.0.0.0
