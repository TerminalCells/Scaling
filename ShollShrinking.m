% ShollShrinking.m

% By Lena Barrett
% lbarrett@princeton.edu
% March 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

clc % Erase command window
clear all % Erase workspace
close all % Close any open figure windows

middleString = 'REPLACE WITH YOUR NICKNAME FOR A GROUP OF TERMINAL CELLS'; % For naming .mat files. I.e., ‘WT’ for the WT group 
graphFileFolder = 'REPLACE WITH FILE PATH TO WHERE EDGE LISTS ARE STORED'; % Path for edge lists
string = 'REPLACE WITH FILE PATH TO WHERE SHRUNKEN EDGE LISTS WILL BE STORED';

addpath(graphFileFolder); % Add graphFileFolder path for MATLAB to search

d = dir(graphFileFolder); % Retrieves structure object listing file info
numFile = length(d); % Gets # of files+2

sampleInfo = zeros(numFile-2,3); % Storing information about samples
% col1 = Sample ID
% col2 = Tracheal Metamere
% col3 = Left (0) or Right (1)

for i = 1:numFile % For every file. Make sure to skip hidden files!
    fileName = d(i).name; % Get file name
    sampleInfo(i-2,1) = str2double(extractBefore(fileName,'_Tr')); % Samp ID
end % End of for loop

replicateNum = max(sampleInfo(:,1)); % Obtains highest replicate #
shrink = 350/560; % The "shrink" factor. Desired max radius/Curr max radius

fileNamePlot = zeros(1,replicateNum); % Initializing for file name
for i = 1:numFile
    fileName = d(i).name; % Extract file name
    fileNamePlotTemp = str2double(extractBefore(fileName,'_Tr')); % ID

    nodes = readmatrix(fileName); % Load data from edge list
    
    child = nodes(2:end,1); % Node IDs
    parent = nodes(2:end,2); % Parent IDs
    [rows,cols] = size(nodes); % # of nodes/data fields

    G = graph(parent,child); % Create graph object
    G.Nodes.x = nodes(:,3); % Add x-coordinates to node
    G.Nodes.y = nodes(:,4); % Add y-coordinates
    G.Nodes.z = nodes(:,5); % Add z-coordinates
    try % If each node has a fitted radius
        G.Nodes.r = nodes(:,6); % Add radii to G
    catch % If each node lacks a fitted radius
    end % End of try-catch statement
        
    afterString = extractAfter(fileName,'_Tr'); % Get metamere #
    sampleInfo(i-2,2) = str2double(afterString(1)); % Save metamere

    if afterString(2) == 'L' % If we're looking at a left cell
        j = 1; % Variable for indicating left vs. right terminal cell
    else % If we're looking at a right cell
        j = 2;
    end % End of if statement
    
    index = i-2; % Need to down-increment due to hidden files on com
    row = sampleInfo(index,2)-1; % Row to save data in, in cell arrays
    sampleInfo(index,3) = j-1; % Recording it's right/left cell

    root = [G.Nodes.x(1) G.Nodes.y(1) G.Nodes.z(1)]; % Stalk-terminal junct
    nodeNum = numnodes(G); % # of SNT nodes in the terminal cell
    newCoords = zeros(nodeNum,3); % Matrix 4 holding new coords
    newCoords(1,:) = root; % Keep same root coordinates as old cell
    for k = 2:nodeNum % For all SNT nodes except the root
        point = [G.Nodes.x(k) G.Nodes.y(k) G.Nodes.z(k)]; % Old coord point
        fullDist =  norm(root-point); % Dist from old SNT node to root
        partDist = shrink*fullDist; % Desired dist from old SNT node
        xyd = sqrt((point(1)-root(1))^2+(point(2)-root(2))^2+(point(3)-root(3))^2);
        for l = 1:3 % For each coordinate
            newCoords(k,l) = root(l) + partDist * (point(l)-root(l))/xyd;
        end % End of for loop
    end % End of for loop

    GSave = G; % Copy characteristics of G into a new graph variable
    GSave.Nodes.x = newCoords(:,1); % Save new, "shrunken" x-coordinate
    GSave.Nodes.y = newCoords(:,2); % Save new, "shrunken" y-coordinate
    GSave.Nodes.z = newCoords(:,3); % Save new, "shrunken" z-coordinate
    save(strcat(string,erase(fileName,'.txt'),'.mat'),'GSave') % Save .mat file; hull area
end % End of for loop

%%
% For visualizing the coordinates of newly shrunken terminal cell, if
% desired, with coordinates of old terminal cell.

figure % Opens a new figure window
    plot3(G.Nodes.x,G.Nodes.y,G.Nodes.z,'o') % Plot each SNT node as circle
    hold on % To ensure that next plot is put on the same figure
    plot3(newCoords(:,1),newCoords(:,2),newCoords(:,3),'ro') % New coords
