// DeleteItAll.ijm

// By Lena Barrett
// lbarrett@princeton.edu
// April 2022
// Shvartsman Laboratory
// Princeton University

// Copyright 2022 Lena Barrett, Princeton University

// The custom IJ1 Macro script to be used after running AllROIs.ijm. First, you open another channel of the image stack you ran AllROIs.ijm on. For example, this could be the nuclear or UV signal of the same image. Then, you run this script, which will remove almost all areas outside of that slice's ROI. You will need to delete the very first slice after running this program.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

for (n=1; n<=nSlices; n++) // For almost every slice in the active, selected z-stack
{
	setSlice(n); // Make the slice selected/active
	roiManager("Select", n); // Selects the ROI for the active slice from the ROI Manager
	run("Clear", "slice"); // Sets all pixels outside of the ROI to have an intensity value of 0
} // End of for loop
