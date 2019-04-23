# Logic on how to launch container depending on a few pre-defined image variables

#Check if Arduino enabled, if yes upload new arduino control
#diff /root/Euplotid/arduino_code/euplouino.ino /data/euplouino.ino || PROGRAMMER=1
if [ "${EUPLO_ARDUINO:-}" == "true" ] ; then
  cd /root/Euplotid/arduino_code/
  make upload && cp euplouino.ino /data/
  unset PROGRAMMER
fi

#If camera installed fire up the RPi cam interface
if [ "${EUPLO_RPICAM:-}" == "true" ]; then
	cd /usr/src/app/RPi_Cam_Web_Interface && ./start.sh
fi

#Set up captive portal for wifi if installed and needed
if [ "${EUPLO_WIFI:-}" == "true" ] ; then
	export DBUS_SYSTEM_BUS_ADDRESS=unix:path=/host/run/dbus/system_bus_socket
	sleep 1
	iwgetid -r

	if [ $? -eq 0 ] ; then
    	printf 'Skipping WiFi Connect\n'
	else
		printf 'Starting WiFi Connect\n'
		cd /usr/src/app/ && ./wifi-connect --clear=false
    fi
fi

# Deploy euplotid app, few scenarios depending on:
# OS: [debian, alpine, ubuntu]
# Arch: [armv6, armv7, aarch64, x86]
# Deploy: [true, false]

if [ "${EUPLO_DEPLOY:-}" == "true" ] ; then 
	python3 /app/applotid.py &
fi

#Dont deploy just run terminal for nanotid
if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" == "nanotid" ] ; then
	if [ "${EUPLO_OS:-}" == "alpine" ] ; then
		/bin/ash
	fi
	if [ "${EUPLO_OS:-}" == "debian" ] ; then
		/bin/bash
	fi
fi
#Run jupyter as IDE for bigger images
if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" == "minitid" ] ; then
	jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0

fi
if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" == "eulertid" ] ; then
	jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0

fi
if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" == "euplotid" ] ; then
	jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0

fi
if [ "${EUPLO_DEPLOY:-}" == "false" ] && [ "${EUPLO_IMAGE:-}" == "megatid" ] ; then
	jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0

fi

