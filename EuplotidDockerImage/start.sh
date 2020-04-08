#!/bin/bash
#Activate default environment
source activate base

# Deploy euplotid app, few scenarios depending on:
# OS: [debian, alpine, ubuntu]
# Arch: [armv6, armv7, aarch64, x86]
# Deploy: [true, false]

#If deploying w/ frontend or not
if [ "${EUPLO_DEPLOY:-}" == "true" ] ; then
	python3 /app/applotid.py
fi

#Run jupyter as IDE
if [ "${EUPLO_DEPLOY:-}" == "false" ] ; then
	jupyter lab --allow-root --port=$JUPYTER_PORT --no-browser --ip=0.0.0.0
fi
