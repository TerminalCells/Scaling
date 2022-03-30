# ExtractingGraph.py

# By Lena Barrett
# lbarrett@princeton.edu
# March 2022
# Shvartsman Laboratory
# Princeton University

# Copyright 2022 Lena Barrett, Princeton University

# Documentation Resources: https://imagej.net/SNT:_Scripting
# Latest SNT API: https://morphonets.github.io/SNT/

# The custom Jython script used to convert traces stored in SNT’s .traces format into edge lists stored in .txt files.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os # For querying directories on your machine
from ij import IJ, ImagePlus # For dealing with GFP micrographs
from sc.fiji.snt import Path, Tree, SNTService, PathFitter, SNT # SNT-specific packages

DIR = ‘PATH FOR .TRACES FILES’

traces = os.listdir(DIR) # Store file names
	
for i in traces: # For each .traces file; Be sure to account for any hidden files which may be in the directory!
	tracesFile = DIR + i # Getting file location
	tree = Tree(tracesFile) # Loading .traces file

	GFP = i.replace('traces','tif') # Save GFP file name
	img_path = ‘PATH FOR GFP MICROGRAPHS’ + GFP # Image path
	IJ.open(img_path) # Open image in Fiji
	img = IJ.getImage() # Shouldn't be necessary in theory but I needed it to get img object defined

	pathList = tree.list() # Retrieve all paths from tree

	for path in pathList: # Fitting each path
		fitter = PathFitter(img, path) # Setup PathFitter  object for GFP micrograph and the path
		fitter.setScope(PathFitter.RADII_AND_MIDPOINTS) # Fit radii AND also optimize node position to be more in the midpoint
		fitter.setMaxRadius(10) # Do not fit if structures have a radius >10 um
		fitter.setReplaceNodes(True) # Replace path's nodes with the “fitted” nodes
		fitter.call() # Perform the fitting

	img.close() # Close open image
	
	os.chdir('PATH FOR EDGE LISTS’) # Directory for saving edge
	
	txt = i.replace('traces','txt') # Setup .txt file to store edge lists

	allNodes = tree.getNodesAsSWCPoints() # Save nodes as SWCPoint objects

	save_file = open(txt, 'w') # Create empty .txt file
	
	for t in allNodes: # For each node
		save_file.write(str(t.id) + ' ' + str(t.parent) + ' ' + str(t.getX()) + ' ' + str(t.getY()) + ' ' + str(t.getZ()) + ' ' + str(t.radius) + '\n') # Write node ID, parent ID, X/Y/Z coordinates...
	
	save_file.close() # Stop writing .txt file
