#!/usr/bin/python
import sys
import logging
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, "/aloy/data/web/sbnb_web/sbnb_web-7.59/covid19/")
from web.app import app as application
