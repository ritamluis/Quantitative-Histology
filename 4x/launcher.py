# python launcher.py N
# Spawns N processes of batch.py, dividing images in ./data into N batches
# Processes are run in background

import os
import sys
import pickle
from numpy import ceil



# Remove any old pickle files
L = os.listdir(".")
for l in L :
	if l.endswith(".p") :
		os.remove(l)


# Create results directory if it doesn't exist
if not os.path.exists("results/") :
	os.makedirs("results/")


# Read files
files = []
directory = "data/"
for f in os.listdir(directory) :
    if f.endswith(".jpg") : # if it's a jpg
        if not f.endswith("_processed.jpg") : # and isn't a processed image
            files.append(directory + f) # then add it to the list to be processed


# Confirm files
print "%d files, requires %d simultaneous runs." % (len(files), ceil(float(len(files)) / float(sys.argv[1])))



# Generate a list of files for each batch.py to process, and spawn processes
for i in range(int(sys.argv[1])) :
	F = open("filelist%i.p" % i, "w")
	l = [f for f in files[i :: int(sys.argv[1])]]
	pickle.dump(l, F)
	F.close()
	#os.system("python batch.py %i &" % i)
