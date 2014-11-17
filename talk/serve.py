import os
import sys
import glob
import numpy as np
from subprocess import call
from IPython.nbconvert.postprocessors import ServePostProcessor

# Grab relevant file
file = glob.glob("%s*.ipynb" % sys.argv[1])[0]

# Generate Reveal.js slides
call(["ipython", "nbconvert", file, "--to", "slides"])

# Serve
server = ServePostProcessor(port = 8000 + np.random.randint(low=10, high=1000))
server(file.split(".ipynb")[0] + ".slides.html?theme=serif")
