% ConvexHull_Graph.m

% By Lena Barrett
% lbarrett@princeton.edu
% March 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

% This script performs more complex calculations to prepare for scaling
% analyses

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clc % Erase command window
clear all % Erase workspace
close all % Close any open figure windows

calibration = REPLACE WITH YOUR MICROSCOPES CALIBRATION; % um/pixels

middleString = 'REPLACE WITH GROUP NICKNAME FOR FILES'; % For naming .mat files. I.e., ‘WT’ for the WT group
graphFileFolder = 'REPLACE WITH DIRECTORY YOU WANT TO SAVE EDGE LISTS IN’; % Path for edge lists

addpath(graphFileFolder); % Add graphFileFolder path for MATLAB to search

d = dir(graphFileFolder); % Retrieves structure object listing file info
numFile = length(d); % Gets # of files+2

sampleInfo = zeros(numFile-2,3); % Storing information about samples
% col1 = Sample ID
% col2 = Tracheal Metamere
% col3 = Left (0) or Right (1)

for i = 3:numFile % For every file, skipping the first two which are N/A
    fileName = d(i).name; % Get file name
    sampleInfo(i-2,1) = str2double(extractBefore(fileName,'_Tr')); % Samp ID
end % End of for loop

replicateNum = max(sampleInfo(:,1)); % Obtains highest replicate #

aPlot{1} = avPlot{1}; % Initializes area of convex hull matrix for left C's
aPlot{2} = avPlot{2}; % Initializes area of convex hull matrix for right Cs
LPlot{1} = aPlot{1}; % Initializes sum of 2D branch lengths for left TCs
LPlot{2} = aPlot{1}; % Initializes sum of 2D branch lengths for right TCs
DistPlot{1} = LPlot{1}; % Initializes projection error matrix for left TC's
DistPlot{2} = DistPlot{1}; % Initializes projection error for right TC's
distPlot = cell(1,2); % Initializes cell array for average void radii
fileNamePlot = zeros(1,replicateNum); % For organizing file names by reps

for i = 1:numFile % For every file in directory; be careful about hidden files!!
    fileName = d(i).name; % Extract file name

    nodes = readmatrix(fileName); % Load data from edge list
    
    child = nodes(2:end,1); % Node IDs
    parent = nodes(2:end,2); % Parent IDs
    [rows,cols] = size(nodes); % # of nodes/data fields

    G = graph(parent,child); % Create graph object
    G.Nodes.x = nodes(:,3); % Add x-coordinates to node
    G.Nodes.y = nodes(:,4); % Add y-coordinates
    G.Nodes.z = nodes(:,5); % Add z-coordinates
    try % If radii have been fit to this terminal cell
        G.Nodes.r = nodes(:,6); % Save fitted radius for each node
    catch % If radii have not been fit to this terminal cell
    end % End of try-catch statement
        
    afterString = extractAfter(fileName,'_Tr'); % Get tracheal metamere #
    sampleInfo(i-2,2) = str2double(afterString(1)); % Save trach. metamere

    if afterString(2) == 'L' % If we're looking at a left cell
        j = 1; % Variable for indicating left vs. right terminal cell
    else % If we're looking at a right cell
        j = 2;
    end % End of if statement
    
    index = i-2; % Need to down-increment by two due to hidden files in my computer
    row = sampleInfo(index,2)-1; % Row to save data in, in cell arrays
    sampleInfo(index,3) = j-1; % Recording whether it's a right/left cell

    depth = max(G.Nodes.z) - min(G.Nodes.z); % Spread of cell along z axis

    xyz = [G.Nodes.x G.Nodes.y G.Nodes.z]; % Combining coords into matrix
    xyzc = mean(xyz); % Saving mean x, y, and z coords
    xyztrans = xyz - xyzc; % Translating points about centroid
    [~,~,V] = svd(xyztrans,0); % Singular value decomp
        % SVD: Factorization (eigendecomposition: representing in terms of
        % eigenvals/vecs)
        % V is an orthogonal matrix; full of orthogonal and normalized vecs
        % Therefore, V is a vector basis
    xyz2d = xyz; % Initializes coordinates for the flattened terminal cell
    numCoords = length(xyz); % Saves # of SNT nodes
    Dist = zeros(numCoords,1); % Initializes Dist matrix with zeros
    normal = [V(1,3) V(2,3) V(3,3)]; % Vector normal to plane of flattened cell
    intercept = -xyzc*normal'; % For "intercept" in plane equation
    [xx,yy] = meshgrid(1:10:160,1:10:140); % Defines coords covered by plane
    z = (-normal(1)*xx - normal(2)*yy - intercept)/normal(3); % Solve for z in plane eq
    xyz(:,3) = xyz(:,3) - z(1,1); % Transforming z-coordinates 
    for k = 1:numCoords % For every SNT node
        baseD = dot(xyz(k,:),V(:,3)); % Length along plane normal to plane
        Dist(k) = abs(baseD); % Abs value of length along plane norm 2 plane
        xyz2d(k,:) = xyz(k,:) - (baseD*V(:,3))'; % Save new coordinate
    end % End of for loop

    DistNorm = Dist/(max(xyz(:,1))-min(xyz(:,1))); % Normalization by LR sp

    G2d = G; % Initializing graph variable containing the flattened cell
        G2d.Nodes.x = xyz2d(:,1); % Save x coordinates
        G2d.Nodes.y = xyz2d(:,2); % Save y coordinates
        G2d.Nodes.z = xyz2d(:,3); % Save z coordinates

    D = degree(G2d); % Degree of each node
    DLength = length(D); % # of nodes
    numEdges = numedges(G2d); % # of edges/branches in the terminal cell
    for k = 1:numEdges % For every edge/branch
        X = [G2d.Nodes.x(G2d.Edges.EndNodes(k,1)),G2d.Nodes.y(G2d.Edges.EndNodes(k,1)),G2d.Nodes.z(G2d.Edges.EndNodes(k,1))
             G2d.Nodes.x(G2d.Edges.EndNodes(k,2)),G2d.Nodes.y(G2d.Edges.EndNodes(k,2)),G2d.Nodes.z(G2d.Edges.EndNodes(k,2))];
        G2d.Edges.Weight(k) = pdist(X); % Save length between SNT nodes
        X = [G.Nodes.x(G.Edges.EndNodes(k,1)),G.Nodes.y(G.Edges.EndNodes(k,1)),G.Nodes.z(G.Edges.EndNodes(k,1))
             G.Nodes.x(G.Edges.EndNodes(k,2)),G.Nodes.y(G.Edges.EndNodes(k,2)),G.Nodes.z(G.Edges.EndNodes(k,2))];
        G.Edges.Weight(k) = pdist(X); % Save length between 3D SNT nodes
    end % End of for loop
    
    [k,aPlot{j}(row,sampleInfo(index,1))] = convhull(xyz2d(:,1),xyz2d(:,2));

    fig1 = figure('visible','off'); % Setting up figure to hold hidden fig
        fill(xyz2d(k,2),xyz2d(k,1),'r'); % Drawing 2D convex hull perimeter
        set(gca,'DataAspectRatio',[1 1 1]) % Set uniform data aspect ratio
        xlim([0 1000]) % Set x axis limits from 0 to 1000 um
        ylim([0 1000]) % Set y axis limits from 0 to 1000 um
        axis off; % Hide both axes in figure
        frm = getframe(gcf); % Save figure
        img = frame2im(frm); % Convert figure to image
        img = rgb2gray(img); % Convert color image to grayscale image
        imgMask = imbinarize(img); % Binarize image
        imgMask = imcomplement(imgMask); % Invert binarized convex hull
    fig2 = figure('visible','off'); % Setting up figure for holding edges
        set(gca,'DataAspectRatio',[1 1 1])
        xlim([0 1000])
        ylim([0 1000])
        try % If radii are known for SNT nodes
            for l = 1:numEdges % For every edge
                r = mean([nodes(G.Edges.EndNodes(l,1),6) nodes(G.Edges.EndNodes(l,2),6)]);
                line([xyz2d(G.Edges.EndNodes(l,1),2) xyz2d(G.Edges.EndNodes(l,2),2)],[xyz2d(G.Edges.EndNodes(l,1),1) xyz2d(G.Edges.EndNodes(l,2),1)],'LineWidth',r)
            end % End of for loop
        catch % If radii are not known for SNT nodes
            for l = 1:numEdges
                line([xyz2d(G.Edges.EndNodes(l,1),2) xyz2d(G.Edges.EndNodes(l,2),2)],[xyz2d(G.Edges.EndNodes(l,1),1) xyz2d(G.Edges.EndNodes(l,2),1)],'LineWidth',1)
            end % End of for loop
        end % End of try-catch statement
        axis off;
        frm = getframe(gcf);
        img = frame2im(frm);
        img = rgb2gray(img);
        imgBranches = imbinarize(img); 
        imgBranches = imcomplement(imgBranches);
        imshow(imgBranches) % Show binarized image of branches on figure
        D = bwdistgeodesic(imgMask,imgBranches); % Distance transform of im
        DClean = D(~isnan(D)); % Remove NaN values, where no distance trans
        DClean = nonzeros(DClean); % Remove 0 values, where branches were

    LPlot{j}(row,sampleInfo(index,1)) = sum(G2d.Edges.Weight); % 2D Length
    DistPlot{j}(row,sampleInfo(index,1)) = mean(DistNorm); % Average proj E
    distPlot{j}(row,sampleInfo(index,1)) = mean(DClean)*calibration; % Rv
end % End of for loop
    
aPlotFlipped = aPlot; % Initializes array with rearranged rows
aPlotPlot = cell(1,2); % Initializes hull area array for scatter plots
LPlotPlot = cell(1,2); % Initializes 2D lengths array for scatter plots
LPlot3DPlot = cell(1,2); % Initializes 2D lengths array; scatter plots
DistPlotFlipped = DistPlot; % Initializes cell array for rearranged rows
DistPlotPlot = cell(1,2); % Initializes cell array for fully rearranged mat
distPlotFlipped = distPlot; % Initializes cell array for rearranged rows
distPlotPlot = cell(1,2); % Initializes cell array for full arranged R .mat
for j = 1:2 % For both left and right cells
    for k = 1:8 % For each tracheal metamere
        aPlotFlipped{j}(k,:) = aPlot{j}(9-k,:); % Rearranges/flips rows
        LPlotFlipped{j}(k,:) = LPlot{j}(9-k,:); % Rearranges/flips rows
        DistPlotFlipped{j}(k,:) = DistPlot{j}(9-k,:); % Rearrange rows
        distPlotFlipped{j}(k,:) = distPlot{j}(9-k,:); % Rearrange rows
    end % End of for loop
    aPlotPlot{j} = aPlotFlipped{j}(:,any(aPlotFlipped{j}))/1e6; % mm
    LPlotPlot{j} = LPlotFlipped{j}(:,any(LPlotFlipped{j}))/1e3; % mm
    DistPlotPlot{j} = DistPlotFlipped{j}(:,any(DistPlotFlipped{j})); % %
    distPlotPlot{j} = distPlotFlipped{j}(:,any(distPlotFlipped{j})); % %
end % End of for loop

string = 'REPLACE WITH DIRECTORY WHERE YOU’LL SAVE .MAT VARIABLES'; % Base path
bigString = strcat(string,middleString); % Base path + sample type descrip
save(strcat(bigString,'_LvV.mat'),'avPlotPlot') % Save .mat file; hull vol
save(strcat(bigString,'_L2d.mat'),'LPlotPlot') % Save .mat file; 2D lengths
save(strcat(bigString,'_Dist.mat'),'DistPlotPlot') % Save .mat file; Distan
save(strcat(bigString,'_DtoBranch.mat'),'distPlotPlot') % Save .mat file
