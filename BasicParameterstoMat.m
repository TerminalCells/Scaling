% BasicParameterstoMat.m

% By Lena Barrett
% lbarrett@princeton.edu
% March 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

% This script performs simple calculations using terminal cells traced in
% SNT, processed using the ExtractingGraph.py script, then processed using
% the DendrogramFromFile.m script. Outputs are .mat variables containing
% ordered (by left/right cell (cell column), by tracheal metamere
% (matrix rows), and by replicate ID (matrix columns)) values for the
% number of endpoints, nodes, and edges AND for the sum of edge lengths
% grouped by Horton-Strahler Order 1 (first row in cell array) and Orders
% >1 (second row in cell array).

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clc % Erase command window
clear all % Erase workspace
close all % Close any open figure windows

middleString = 'INSERT NICKNAME FOR GROUP OF TERMINAL CELLS'; % For naming .mat files
graphFileFolder = 'INSERT PATH TO RETRIEVE .TXT EDGE LISTS'; % Path for edge lists
matFileFolder = 'INSERT PATH TO RETRIEVE .MAT FILES CONTAINING INFORMATION ABOUT TERMINAL CELLS'; % Path for .mat

addpath(graphFileFolder); % Add graphFileFolder path for MATLAB to search
addpath(matFileFolder);

d = dir(graphFileFolder); % Retrieves structure object listing file info
numFile = length(d); % Gets # of files+2

sampleInfo = zeros(numFile-2,3); % Storing information about samples
% col1 = Sample ID
% col2 = Tracheal Metamere
% col3 = Left (0) or Right (1)

for i = 3:numFile % For every file, skipping the first two which are N/A
    fileName = d(i).name; % Get file name
    try % Should work if there are no hidden files
        sampleInfo(i-2,1) = str2num(extractBefore(fileName,'_Tr')); % Samp ID
    catch % Skip if this is a hidden file
    end % End of try-catch statement
end % End of for loop

replicateNum = max(sampleInfo(:,1)); % Obtains highest replicate #

endptNumPlot = cell(1,2); % Initializes cell array for holding endpoint #
nodePlot = cell(1,2); % Initializes cell array for holding node # by L/R Ce
edgePlot = cell(1,2); % Initializes cell array for holding edge # by L/R Ce
HSPlot = cell(2,2); % Initializes cell array 4 holding sum of edge L by HS
    for l = 1:2 % For left and right cells
        for k = 1:2 % For Horton-Strahler Order 1 and >1
            HSPlot{l,k} = zeros(8,replicateNum); % Initializes matric w/ 0s
        end % End of for loop
    end % End of for loop
for i = 1:numFile % Be sure to account for hidden files!!
    fileName = d(i).name; % Extract file name
    old = load(strcat(matFileFolder,strcat(erase(fileName,'txt'),'mat'))); % Loading .mat file into workspace
    D = degree(old.G); % Getting degree of each node in graph

    afterString = extractAfter(fileName,'_Tr'); % Removing part of string
    sampleInfo(i-2,2) = str2num(afterString(1)); % Saving tracheal metamere
    if afterString(2) == 'L' % If we're looking at a left cell
        j = 1; % Variable for indicating left vs. right terminal cell
    else % If we're looking at a right cell
        j = 2;
    end % End of if statement

    index = i-2; % Need to down-increment by two due to non-files
    row = sampleInfo(index,2)-1; % Row to save data in, in cell arrays
    sampleInfo(index,3) = j-1; % Recording whether it's a right/left cell

    edgeNum = numedges(old.G); % Number of edges
    for k = 1:edgeNum % For every edge
        if old.G.Edges.HSOrder(k) == 1 % For edge of Horton-Strahler Ord 1
            HSPlot{j,1}(row,sampleInfo(index,1)) = HSPlot{j,1}(row,sampleInfo(index,1)) + old.G.Edges.Weight(k);
        else % For edges of Horton-Strahler Order >1
            HSPlot{j,2}(row,sampleInfo(index,1)) = HSPlot{j,2}(row,sampleInfo(index,1)) + old.G.Edges.Weight(k);
        end % End of if statement
    end % End of for loop

    endptNumPlot{j}(row,sampleInfo(index,1)) = sum(D == 1) - 1; % Calc Ne

    nodePlot{j}(row,sampleInfo(index,1)) = numnodes(old.G) - 1; %Calc node#

    edgePlot{j}(row,sampleInfo(index,1)) = numedges(old.G); % Calcul edge #
end % end of for loop

endptNumPlotFlipped = endptNumPlot; % For reorganizing rows of matrices
endptNumPlotPlot = cell(1,2); % Initializing fully organized matrices
nodePlotFlipped = nodePlot; % For reorganizing rows of matrices
nodePlotPlot = cell(1,2); % Initializing fully organized matrices
edgePlotFlipped = edgePlot; % For reorganizing rows of matrices
edgePlotPlot = cell(1,2); % Initializing fully organized matrices
HSPlotFlipped = HSPlot; % For reorganizing rows of matrices
HSPlotPlot = cell(2,2); % Initializing fully organized matrices
for j = 1:2 % For left and right cells
    for k = 1:8 % For each tracheal metamere (body segment)
        endptNumPlotFlipped{j}(k,:) = endptNumPlot{j}(9-k,:); % "Flip" rows
        nodePlotFlipped{j}(k,:) = nodePlot{j}(9-k,:); % Flip/reorder rows
        edgePlotFlipped{j}(k,:) = edgePlot{j}(9-k,:); % Flip/reorder rows
    end % End of for loop
    endptNumPlotPlot{j} = endptNumPlotFlipped{j}(:,any(endptNumPlotFlipped{j}));
    nodePlotPlot{j} = nodePlotFlipped{j}(:,any(nodePlotFlipped{j})); % Node
    edgePlotPlot{j} = edgePlotFlipped{j}(:,any(edgePlotFlipped{j})); % Edge
end % End of for loop
for l = 1:2 % For Horton-Strahler Order 1 and >1
    for j = 1:2 % left and right cells
        for k = 1:8 % for each tracheal metamere (body segment)
            HSPlotFlipped{j,l}(k,:) = HSPlot{j,l}(9-k,:); % Reorder rows
        end % End of for loop
        HSPlotPlot{j,l} = HSPlotFlipped{j,l}(:,any(HSPlotFlipped{j,l}));
    end % End of for loop
end % End of for loop

string = 'Y:\Lena\Analyses\ExtraAnalyses\.mat Variables\'; % Base path
bigString = strcat(string,middleString); % Base path + sample type descrip
save(strcat(bigString,'_EndPt.mat'),'endptNumPlotPlot') % Save .mat file
save(strcat(bigString,'_Nodes.mat'),'nodePlotPlot') % Save .mat file
save(strcat(bigString,'_edges.mat'),'edgePlotPlot') % Save .mat file
save(strcat(bigString,'_HS.mat'),'HSPlotPlot') % Save .mat file
