% ShollCorrect.m

% By Lena Barrett
% lbarrett@princeton.edu
% March 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

% The custom MATLAB script used to conduct Sholl analyses.

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clc % Erase command window
clear all % Erase workspace
close all % Close any open figure windows

m = 0; % If true, use the shrunken edge lists; 1 = true, 0 = false

middleString = 'REPLACE WITH YOUR NICKNAME FOR A GROUP OF TERMINAL CELLS'; % For naming .mat files. I.e., ‘WT’ for the WT group
if m == 0 % For using “normal” edge lists
    graphFileFolder = 'REPLACE WITH DIRECTORY OF YOUR NORMAL EDGE LISTS'; % Directory to retrieve edge lists
else % For using shrunken edge lists
    graphFileFolder = 'REPLACE WITH DIRECTORY OF YOUR SHRUNKEN EDGE LISTS'; % Directory to retrieve shrunken edge lists; should only be relevant for WT late L3 ones
end % End of if statement

addpath(graphFileFolder); % Add graphFileFolder path for MATLAB to search
d = dir(graphFileFolder); % Retrieves structure object listing file info
numFile = length(d); % Gets # of files+2

sampleInfo = zeros(numFile-2,3); % Storing information about samples
% col1 = Sample ID
% col2 = Tracheal Metamere
% col3 = Left (0) or Right (1)

for i = 1:numFile % For every file, BE SURE TO SKIP HIDDEN FILES
    fileName = d(i).name; % Get file name
    sampleInfo(i-2,1) = str2double(extractBefore(fileName,'_Tr')); % Samp ID
end % End of for loop

replicateNum = max(sampleInfo(:,1)); % Obtains highest replicate #
binSize = 10; % Concentric circles are spaced this many microns apart

%%
% After running the top section to initialize, run this section if the .mat
% Sholl file needs to be updated.

sholl = cell(1,2); % Initializes matrices holding each cell's Sholl plot
sholl{1} = cell(8,replicateNum); % 8 = # tracheal metameres, for left cells
sholl{2} = sholl{1}; % Copy of sholl{1} for right cells; for initialization
fileNamePlot = zeros(1,replicateNum); % For keeping track of replicates
for i = 3:numFile % For every file in directory
    fileName = d(i).name; % Extract file name
    fileNamePlotTemp = str2double(extractBefore(fileName,'_Tr')); % ID

    afterString = extractAfter(fileName,'_Tr'); % Get tracheal metamere #
    row = 10 - str2double(afterString(1)); % Save trach. metamere
    if afterString(2) == 'L' % If we're looking at a left cell
        leftoright = 1; % Variable for indicating left v right cell
    else % If we're looking at a right cell
        leftoright = 2;
    end % End of if statement

    try % First, try loading terminal cell network from edge list
        nodes = readmatrix(fileName); % Load data from edge list
        
        child = nodes(2:end,1); % Node IDs
        parent = nodes(2:end,2); % Parent IDs
        [rows,cols] = size(nodes); % # of nodes/data fields

        G = graph(parent,child); % Create graph object
        G.Nodes.x = nodes(:,3); % Add x-coordinates to node
        G.Nodes.y = nodes(:,4); % Add y-coordinates
        G.Nodes.z = nodes(:,5); % Add z-coordinates
    catch % Next, load terminal cell network from .mat graph file
        loaded = load(strcat(graphFileFolder,fileName)); % Load from .mat
        G = loaded.GSave; % Rename graph
    end % End of try-catch statement

    numEdges = numedges(G); % # of edges (between SNT nodes)
    for j = 1:numEdges % For every SNT node edge
        [~,G.Edges.startDist(j)] = dsearchn([G.Nodes.x(1) G.Nodes.y(1) G.Nodes.z(1)],[G.Nodes.x(G.Edges.EndNodes(j,1)) G.Nodes.y(G.Edges.EndNodes(j,1)) G.Nodes.z(G.Edges.EndNodes(j,1))]);
        [~,G.Edges.endDist(j)] = dsearchn([G.Nodes.x(1) G.Nodes.y(1) G.Nodes.z(1)],[G.Nodes.x(G.Edges.EndNodes(j,2)) G.Nodes.y(G.Edges.EndNodes(j,2)) G.Nodes.z(G.Edges.EndNodes(j,2))]);
    end % End of for loop

    bins = 0:binSize:ceil(max(G.Edges.endDist)/binSize)*binSize; % Array
    bins = vertcat(bins,zeros(1,max(bins)/binSize+1)); % Add histogram line
    lengthBins = length(bins(1,:)); % # of bins
    
    for j = 2:lengthBins % For each bin
        batch1 = find(G.Edges.endDist > bins(1,j)); % Edges w endnode < bin
        batch2 = find(G.Edges.startDist < bins(1,j)); % Startnode > bin
        [val,~] = intersect(batch1,batch2); % Shared edges have intersects
        bins(2,j) = length(val); % # of branch intersections for the bin
    end % End of for loop

    sholl{leftoright}{row,fileNamePlotTemp} = bins(2,:); % Save Sholl
end % End of for loop

string = 'Y:\Lena\Analyses\ExtraAnalyses\.mat Variables\'; % Base dir path
bigString = strcat(string,middleString); % Base path + sample type descrip
if m == 0
    save(strcat(bigString,'_Sholl.mat'),'sholl') % Save .mat file; all individual “normal” Sholl analyses
else
    save(strcat(bigString,'_ShollShrunk.mat'),'sholl') % Save .mat file; individual “shrunken” Sholl analyses
end % End of if statement

%%
% Run after running the first section and if the Sholl .mat file is ready
% to combine all individual Sholl profiles into one representative Sholl
% curve

clc % Clears command window

if m == 0
    loaded = load('INSERT WITH PATH TO .MAT FILE CONTAINING ALL INDIVIDUAL NORMAL SHOLL ANALYSES'); % Load .mat variable in MATLAB workspace 
    sholl = loaded.sholl; % Rename Sholl variable
else
    loaded = load('INSERT WITH PATH TO .MAT FILE CONTAINING ALL INDIVIDUAL SHRUNKEN SHOLL ANALYSES');
    sholl = loaded.sholl; % Rename Sholl variable
end % End of if statement

maxBin = cell(1,2); % Initializes for holding # of bins for each cell
maxBinMax = zeros(1,2); % Initializes for holding maximum # of bins 
for i = 1:2 % For left and right cells
    sholl{i} = horzcat(sholl{i},sholl_E203Kctrl{i}); % Combine two WT group
    maxBin{1,i} = cellfun(@length,sholl{i}); % Get length of each matrix
    maxBinMax(1,i) = max(max(maxBin{i})); % Get max # of bins for L/R cells
end % End of for loop

maxBinMaxMax = max(maxBinMax); % Get max # of bins for ALL cells

bins = 0:binSize:maxBinMaxMax*binSize; % Make std array of bins
bins = vertcat(bins,zeros(1,maxBinMaxMax+1)); % Add histogram line

dif = cell(1,2); % Initializes variable for # of "missing" replicates
for i = 1:2 % For left and right cells
    dif{i} = maxBinMaxMax - maxBin{i}; % # bins each indiv prof is missing
end % End of for loop

shollPad = cell(1,2); % Initializes copy of sholl, but with padded zeros
n = cell(1,2); % Initializes variable for holding # of animals in data
stdev = n; % Initializes variable for storing std dev of each sholl distrib
shollDist = cell(1,2); % Initializes variable for mean of each sholl distri
replicateNum = max(max(size(sholl{1}),size(sholl{2}))); % Max replicate val
for i = 1:2 % For left and right cells
    missing = cellfun(@isempty,sholl{i}); % Find places where data missing
    n{i} = zeros(8,1); % Fill with matrices of 0s: rows are metameres
    stdev{i} = cell(8,1); % Fill with cell array: rows are metameres
    shollDist{i} = cell(8,1); % Fill with cell array: rows are metameres 
    for j = 1:8 % For each tracheal metamere
        shollDist{i}{j} = zeros(1,maxBinMaxMax); % Fill with matrix of 0s
        stdev{i}{j} = shollDist{i}{j}; % Fill with matrix of 0s: cols= bins
        for k = 1:replicateNum % For each replicate
            shollPad{i}{j,k} = padarray(sholl{i}{j,k},[0 dif{i}(j,k)],'post');
        end % End of for loop
        n{i}(j,1) = replicateNum - sum(missing(j,:)); % Subtract empties 
        for k = 1:replicateNum % For each replicate
            if isempty(shollPad{i}{j,k}) == 0 % If data are not missing
                shollDist{i}{j} = shollDist{i}{j} + shollPad{i}{j,k}; % Sum
            end % End of if statement
        end % End of for loop
        shollDist{i}{j} = (shollDist{i}{j}./n{i}(j))'; % Sum/# of animals
        for k = 1:replicateNum % For each replicate
            if isempty(shollPad{i}{j,k}) == 0 % If data are not missing
                for l = 1:maxBinMaxMax % For each bin
                    stdev{i}{j}(l) = stdev{i}{j}(l) + (shollPad{i}{j,k}(l) - shollDist{i}{j}(l))^2;
                end % End of for loop
            end % End of if statement
        end % End of for loop
        stdev{i}{j} = sqrt(stdev{i}{j}./n{i}(j))'; % Putting sum in formula
    end % End of for loop
end % End of for loop
