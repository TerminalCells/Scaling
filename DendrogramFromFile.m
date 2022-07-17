% DendrogramFromFile.m

% By Lena Barrett
% lbarrett@princeton.edu (Will be decommissioned in May 2024)
% July 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

% The custom MATLAB script used to process edge lists just extracted from
% ExtractingGraph.py in Fiji to further characterize graph variables
% describing terminal cells and create a connectivity graph for each
% terminal cell.

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clc % Clear command window
clear all % Clear workspace
close all % Close any open figures

edgeFileFolder = 'INSERT FILE PATH TO WHERE YOUR EDGE LISTS ARE’;
matFileBase = 'INSERT FILE PATH OF WHERE .MAT FILES OF SIMPLIFIED TERMINAL CELLS WILL BE';
pngFileBase = 'INSERT FILE PATH OF WHERE CONNECTIVITY GRAPHS OF SIMPLIFIED TERMINAL CELLS WILL BE';
addpath(edgeFileFolder); % Have MATLAB search in folder with edge lists

d = dir(edgeFileFolder); % Save file properties in edgeFileFolder directory
numFile = length(d); % # of files

for j = 1:numFile % Change to process correct files; be sure to exclude hidden files!
    fileName = d(j).name; % String with file name, not a long directory
    nodes = readmatrix(fileName); % Load data from edge list
    
    child = nodes(2:end,1); % Node IDs
    parent = nodes(2:end,2); % Parent IDs
    [rows,cols] = size(nodes); % # of nodes/data fields

    G = graph(parent,child); % Create graph object
    G.Nodes.x = nodes(:,3); % Add x-coordinates to node
    G.Nodes.y = nodes(:,4); % Add y-coordinates
    G.Nodes.z = nodes(:,5); % Add z-coordinates
    
    Empty = cols == 6; % For managing files with no radii 
        if Empty == 0 % If radii are missing
        else % If radii are not missing
            G.Nodes.r = nodes(:,6); % Save radii
        end % End of if statement

    Go = G; % Save "complicated" graph object for later

    Do = degree(G); % Degree of each node
    DLength = length(Do); % # of nodes

    Drid = ones(DLength,1); % Marking which nodes have a degree of 2
    for i = 1:DLength % For every node
        if Do(i) == 2 % Want to remove nodes of degree 2
            Drid(i) = 0; % 0 = has a degree of 2
        end % End of if statement
    end % End of for loop

    nodesKeep = sum(Drid); % # of nodes to keep
    numEdges = numedges(G); % # of edges in unsimplified network from SNT
    for k = 1:numEdges % For every edge
        G.Edges.Weight(k) = norm([G.Nodes.x(G.Edges.EndNodes(k,1)) G.Nodes.y(G.Edges.EndNodes(k,1)) G.Nodes.z(G.Edges.EndNodes(k,1))] ...
            - [G.Nodes.x(G.Edges.EndNodes(k,2)) G.Nodes.y(G.Edges.EndNodes(k,2)) G.Nodes.z(G.Edges.EndNodes(k,2))]);
    end % End of for loop

    while numnodes(G) > nodesKeep % While there are still nodes to remove
        nodeNum = numnodes(G); % # of nodes
        D = degree(G); % Degree of nodes
        numEdges = numedges(G);
        for i = 1:nodeNum % For every node
            if D(i) == 2 % For every node to remove
               N = neighbors(G,i); % Determine neighboring nodes
               neighEdges = find(G.Edges.EndNodes==i); % Find node w/ edge
               for l = 1:2 % For both neighboring edges
                   if neighEdges(l) > numEdges % Readjust
                       neighEdges(l) = neighEdges(l) - numEdges;
                   end % End of if statement
               end % End of for loop
               wt = G.Edges.Weight(neighEdges(1)) + G.Edges.Weight(neighEdges(2));
               G = addedge(G,N(1),N(2),wt); % Add edge between neighbors
               G = rmnode(G,i); % Remove node
               break % Start over
            end % End of if statement
        end % End of for loop
    end % End of while loop
     
    numEdges = numedges(G); % # of edges
    numNodes = numnodes(G); % # of nodes
    Pmax = 0; % Initialize longest unweighted path variable
    lastNode = 2; % Node 2 is always a child of Node 1, so start there
    PStore = zeros([1,2]); % Initialize variable to store longest unwt path
    forks = find(Do==3); % Find node IDs with degree of 3
    forksLength = length(forks); % # of fork pts
    compare = zeros(1,forksLength); %Initializing variable for finding node
    children = zeros(2,5); % Initialize variable for storing next five node
    
    for i = 2:numNodes % For each node
        P = shortestpath(G,2,i,'Method','unweighted'); % Get shortest path
        PLength = length(P); % # of nodes in shortest path
        if PLength > Pmax % If we find a longer shortest path
            Pmax = PLength; % Update
            PStore = P; % Update
            lastNode = i; % Update what the "last" node is
        end % End of if statement
    end % End of for loop
        
    lastNum = zeros(1,numNodes-1-Pmax); % Initialize nodes after sorting
    count = 0; % Initialize count/tracker variable
    nodeNumSimple = numnodes(G); % # of nodes in simple graph

    for m = 1:nodeNumSimple % For every node in simple graph
        if D(m) == 3 % If the node is a fork point
            x1 = G.Nodes.x(m); % Save its x-coord
            y1 = G.Nodes.y(m); % Save its y-coord
            z1 = G.Nodes.z(m); % Save its z-coord
            for k = 1:forksLength % For every...fork? Try to find same node
                compare(k) = Go.Nodes.x(forks(k)) == x1; % Find same node?
            end % End of for loop
            fork = forks(compare==1); % Did we find same node complicated?
            forkLength = length(fork); % # of identical nodes by x-coord
            compare = zeros(1,forkLength); % Reinitialize y/n compare vari
            for k = 1:forkLength % For every potential same node
                compare(k) = Go.Nodes.y(fork(k)) == y1; % Same y-coord?
            end % End of for loop
            fork = fork(compare==1); % Did we find same node by x/y coords?
            forkLength = length(fork);
            if forkLength == 1 % If there are no more comparisons to make
            else % If we need to single out a fork still
                compare = zeros(1,forkLength); % Initialize comparison vari
                for k = 1:forkLength % For every fork to compare
                    compare(k) = Go.Nodes.z(fork(k)) == z1; % Match z-coord
                end % End of for loop
                fork = fork(compare==1); % Identify corresponding coords
            end % End of if statement
            try % If there are 4+ nodes down from child nodes
                children(:,1) = find(nodes(:,2)==fork(1)); % Save children of fork
                for k = 2:5 % Finding next four nodes down from children
                    for l = 1:2 % For both children nodes
                        children(l,k) = find(nodes(:,2)==children(l,k-1));
                    end % End of for loop
                end % End of for loop
                end1 = children(1,5); % Fourth node from first child node
                end2 = children(2,5); % Fourth node form second child node
                b = norm([x1 y1 z1] - [Go.Nodes.x(end1) Go.Nodes.y(end1) Go.Nodes.z(end1)]);
                c = norm([x1 y1 z1] - [Go.Nodes.x(end2) Go.Nodes.y(end2) Go.Nodes.z(end2)]);
                a = norm([Go.Nodes.x(end1) Go.Nodes.y(end1) Go.Nodes.z(end1)] - [Go.Nodes.x(end2) Go.Nodes.y(end2) Go.Nodes.z(end2)]);

                G.Nodes.angle(m) = round(acosd((b^2+c^2-a^2)/2/b/c)); %Calc
            catch % If there are <4 nodes down from child nodes
                G.Nodes.angle(m) = NaN; % Not a #/cannot calc local bif ang
            end % End of try-catch statement
        else % No local bifurcation angle can be calculated
            G.Nodes.angle(m) = NaN; % Node is not a fork, is root/endpt
        end % End of if statement
    end % End of for loop
                
    for i = 2:Pmax % For every node in the longest shortest unwt path
        dif = PStore(i)-PStore(i-1); % What's the difference in node IDs?
        if dif == 0 % If difference, do nothing
        else % If there is no difference
           for k = 1:dif-1 % For sorting nodes such that longest is on left
              count = count+1; % For correctly assigning indices
              lastNum(count) = PStore(i)-k; % For nodes outside of long SP
           end % End of for loop
        end % End of if statement
    end % End of for loop
        
    countw = 0; % Initializing another count variable
    while isempty(find(lastNum==0,1)) == 0 % Continuing to gather all nodes
        count = count+1; % Update count variable
        countw = countw+1; % Update other count variable
        lastNum(count) = PStore(end) + countw; % Gathering nodes not in SP
    end % End of while loop

    order = [1 PStore lastNum]; % New order of nodes by node ID

    G = reordernodes(G,order); % Reorder nodes in graph

    for i = 1:numEdges % For every edge
       D = degree(G); % Calculate degree of each simplified node
       startNode = G.Edges.EndNodes(i,1); % Start node of edge
       endNode = G.Edges.EndNodes(i,2); % End node of edge
       startCoords = [G.Nodes.x(startNode) G.Nodes.y(startNode) G.Nodes.z(startNode)];
       endCoords = [G.Nodes.x(endNode) G.Nodes.y(endNode) G.Nodes.z(endNode)];

       if height(G.Edges) == 1 % If the terminal cell has only one branch
           G.Edges.HSOrder(i) = 0; % Stand-in for N/A
       else % For most terminal cells
           DLength = length(D); % Number of nodes
           endptPlace = zeros(DLength,1); % Initialize vari for endpt posit
           for k = 2:DLength % For every node, except the root coords
               endptPlace(k,1) = double(D(k)==1); % 1 if node is endpoint
               if endptPlace(k,1) == 1 % If node k is an endpoint
                   row = find(G.Edges.EndNodes(:,2)==k); % Find other node
                   if isempty(row) == 1 % If there is no other node at edge
                       row = find(G.Edges.EndNodes(:,1)==k); % Try othe col
                   end % End of if statement
                   G.Edges.HSOrder(row) = 1; % Assign Horton-Str Order of 1
               end % End of if statement
           end % End of for loop

           allDone = 0; % While false, Horton-S ordering still not complete
           while allDone == 0 % Until all edges have Horton-Strahler order
               for k = 2:DLength % For every node except the root coords
                   if D(k) == 3 % For every fork point (i.e. not endpoint)
                       NBs(:,1) = neighbors(G,k); % Find neighbor nodes
                       edges = [findedge(G,k,NBs(1,1)); findedge(G,k,NBs(2,1)); findedge(G,k,NBs(3,1))];
                       HSOrders = [G.Edges.HSOrder(edges(1,1)); G.Edges.HSOrder(edges(2,1)); G.Edges.HSOrder(edges(3,1))];
                       zeroCheck = find(HSOrders==0); % Find edges wo HS or
                       if length(zeroCheck) == 1 % If there are unass edges
                           nonZeroCheck = find(HSOrders~=0); % Find edges
                           sameCheck = G.Edges.HSOrder(edges(nonZeroCheck(1))) == G.Edges.HSOrder(edges(nonZeroCheck(2)));
                           if sameCheck == 1 % If child branch hav same ord
                               G.Edges.HSOrder(edges(zeroCheck)) = G.Edges.HSOrder(edges(nonZeroCheck(1))) + 1;
                           else % If child branches have different orders
                               HSHigh = max(G.Edges.HSOrder(edges(nonZeroCheck(1))),G.Edges.HSOrder(edges(nonZeroCheck(2))));
                               G.Edges.HSOrder(edges(zeroCheck)) = HSHigh;
                           end % End of if statement
                       end % End of if statement
                   end % End of if statement
               end % End of for loop
               allDone = min(G.Edges.HSOrder); % Any unassigned edges?
           end % End of while loop
       end % End of if statement

       for k = [startNode,endNode] % For beginning and ending nodes of edge
           xGo = find(G.Nodes.x(k) == Go.Nodes.x); % Find matching nodes, x
           yGo = find(G.Nodes.y(k) == Go.Nodes.y); % Find matching nodes, y
           zGo = find(G.Nodes.z(k) == Go.Nodes.z); % Find matching nodes, z
           nodeGo = intersect(intersect(xGo,yGo),zGo); % Find comp coords
           nodeGoLength = length(nodeGo); % Was 1 coordinate set found?
           if nodeGoLength ~= 1 % If >1 coordinate set was found
                startCoord = [G.Nodes.x(1) G.Nodes.y(1) G.Nodes.z(1)]; % G
                rG = round(norm(startCoord - [G.Nodes.x(k) G.Nodes.y(k) G.Nodes.z(k)]));
                for l = 1:nodeGoLength % For every potential matching set
                    rGo = round(norm(startCoord - [Go.Nodes.x(nodeGo(l)) Go.Nodes.y(nodeGo(l)) Go.Nodes.z(nodeGo(l))]));
                    if rGo == rG % If a match is found in terms of distance
                        nodeGo = nodeGo(l); % Save the matching set
                        break % Prematurely end for loop
                    end % End of if statement
                end % End of for loop
           end % End of if statement
           if k == startNode % If trying to match start coords of edge
               startNodeGo = nodeGo; % Save the matching starting coord set
           else % If trying to match end coordinates of edge
               endNodeGo = nodeGo; % Save the matching ending coord set
           end % End of if statement
       end % End of for loop

       P = shortestpath(Go,startNodeGo,endNodeGo); % Map path in original
       try % If there are radii fitted to the original SNT terminal cell
           G.Edges.stdR(i) = std(Go.Nodes.r(P)); % Record the standard dev
           G.Edges.r(i) = round(geomean(Go.Nodes.r(P)),1); % Record R mean
           G.Edges.minR(i) = min(Go.Nodes.r(P)); % Record minimum R of edge
           G.Edges.maxR(i) = max(Go.Nodes.r(P)); % Record maximum R of edge
       catch % If there are no radii fitted to original SNT terminal cell
       end % End of try-catch statement
    end % End of for loop
    
    baseFile = erase(fileName,'txt'); % Base of file names
    matFile = strcat(matFileBase,baseFile,'mat'); % Constructing .mat name
    save(matFile,'G') % Save .mat file storing simplified graphs
    pngFile = strcat(pngFileBase,baseFile,'png'); % Construct .png name
    
    windowLength = 2000; % Size of window in one length in pixels
    fig = figure('visible','off','Position',[0 0 windowLength windowLength]);
    
        edgelabelstr = string(round(G.Edges.Weight)); % For labeling edges
        if Empty == 0 % If no radii
        else % If there is radii
            edgelabelstr = edgelabelstr + '|' + string(G.Edges.r); % L|radi
        end % End of if statement

        p = plot(G,'EdgeLabel',cellstr(edgelabelstr),'NodeLabel',G.Nodes.angle,'LineWidth',2);
            p.EdgeFontSize = 5; % Set font size of edge labels to 5
            p.NodeColor = 'red'; % Make nodes red
            p.NodeFontSize = 6; % Set node font size to 6
            xTextStart = 0.05; % Start Horton-Strahler labels 5% of x-direc
            yTextStart = 0.95; % Start Horton-Strahler labels 95% of y-dire
            text(xTextStart,yTextStart,'HS Order Color Key','Units','normalized')
            text(xTextStart,yTextStart-0.1,'1, Purple','Color','#7E2F8E','Units','normalized')
            text(xTextStart,yTextStart-0.08,'2, Blue','Color','blue','Units','normalized')
            text(xTextStart,yTextStart-0.06,'3, Green','Color','g','Units','normalized')
            text(xTextStart,yTextStart-0.04,'4, Orange','Color','#EDB120','Units','normalized')
            text(xTextStart,yTextStart-0.02,'5, Magenta','Color','m','Units','normalized')
           
        for i = 1:height(G.Edges) % For every edge
            if G.Edges.HSOrder(i) == 1 % For every edge with HS Order of 1
                highlight(p,G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2),'EdgeColor','#7E2F8E')
            elseif G.Edges.HSOrder(i) == 2 % For every edge w/ HS Order 2
                highlight(p,G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2),'EdgeColor','blue')
            elseif G.Edges.HSOrder(i) == 3
                highlight(p,G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2),'EdgeColor','g')
            elseif G.Edges.HSOrder(i) == 4
                highlight(p,G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2),'EdgeColor','#EDB120')
            else % For every edge with Horton-Strahler Order of 5
                highlight(p,G.Edges.EndNodes(i,1),G.Edges.EndNodes(i,2),'EdgeColor','m')
            end % End of if statement
        end % End of for loop

    saveas(fig,pngFile) % Save figure as .png file as per pngFile direc
end % End of for loop
