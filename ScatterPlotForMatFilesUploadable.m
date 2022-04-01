% ScatterPlotForMatFilesUploadable.m

% By Lena Barrett
% lbarrett@princeton.edu
% April 2022
% Shvartsman Lab
% Princeton University

% Copyright 2022 Lena Barrett, Princeton University

% The custom MATLAB script used to 1. load .mat variables containing
% parameters to use for scaling analyses, 2. Construct scatter plots
% showcasing scaling analyses, and 3. apply linear regression techniques.
% YOU WILL NEED TO ALSO DOWNLOAD the MATLAB function maregress.m from
% https://www.mathworks.com/matlabcentral/fileexchange/27916-maregress.

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

clc % Clear command window
clear all % Clear workspace
close all % Close any open figures

savePath = 'INSERT PATH TO DIRECTORY WHERE YOU WANT TO SAVE .PNG IMAGES OF PLOTS';
addpath('INSERT PATH TO DIRECTORY WHERE THE maregress.m FILE IS LOCATED')
alpha = 0.05; % For the maregress.m function; states p-value significance

% Repeat this code for as many .mat variables you want to load
loaded = load('Y:\Lena\Analyses\ExtraAnalyses\.mat Variables\WT_edges.mat');
variableName = loaded.variableNameInMatFile; % Redefine name of variable
isNZ_variableName{1} = (~variableName{1}==0); % Replaces 0's w/ NaN, left C
isNZ_variableName{2} = (~variableName{2}==0); % Replaces 0's w/ NaN, right

%%
% Run this section after running the first section to load your variables

close all % Close any open figures
clc % Clear command window

position = 1; % Which metamere? 1 = Tr9, 2 = Tr8, 3 = Tr7, 4 = Tr6, etc.
fig9a = figure; % Define figure variable to hold your figure

    rowX = xVari{1}(position,:); % Capturing Tr9 data of left cells for x
    rowY = yVari{1}(position,:); % Capturing Tr9 data of left cells for y
    scatter(rowX(isNZ_x{1}(position,:)),rowY(isNZ_y{1}(position,:)),'filled','MarkerFaceColor','black','MarkerFaceAlpha',0.5);
    set(gca,'xscale','log') % Comment out for non-loglog plot
    set(gca,'yscale','log') % Comment out for non-loglog plot
    ylim([1e-2 10]) % Define limits of y-axis: first number is lower bound
    xlabel('Convex Hull Area (mm^2)') % Remove for log-log scale
    ylabel('Sum of Edge Lengths (mm)'); % Remove for log-log scale
    H = gca; % Capture the current axes into a variable
    H.LineWidth = 1; % Increase line width of x/y axes (next to numbering)
    set(H,'FontSize',12) % We recommend 18 for loglog plots, 12 for linear
    hold on % Need to ensure next "scatter" function does not overwrite
    rowX2 = xVari{2}(1,:); % Capturing Tr9 data of right cells for x
    rowY2 = yVari{2}(1,:); % Capturing Tr9 data of right cells for y
    rowX2(rowX2(isNZ_x{2}(position,:)),rowY2(isNZ_y{2}(position,:)),'filled','MarkerFaceColor','black','MarkerFaceAlpha',0.5);
   
    bigX = log10(horzcat(rowX(isNZ_x{1}(position,:)),rowX2(isNZ_x{2}(position,:))));
    bigY = log10(horzcat(rowY(isNZ_y{1}(position,:)),rowY2(isNZ_y{2}(position,:))));

    mdlL3 = fitlm(bigY,bigX); % Ordinary least squares model regression
    x = linspace(1e-2,0.2); % Vector containing values of x for plotting

    [bL3,bintL3,~,~] = maregress(bigY,bigX,alpha); % Major axis regression
    y = 10^(bL3(1))*x.^(bL3(2)); % Vector containing values of y for MA
    plot(x,y,'k') % Plot MA regression with a solid black line
    yConf = 10^(bintL3(1,1))*x.^(bintL3(2,1)); % y values for confidence In
    plot(x,yConf,':k') % Plot one confidence interval line with dotted blac
    yConf = 10^(bintL3(1,2))*x.^(bintL3(2,2)); % y values for other confide 
    plot(x,yConf,':k') % Plot other confidence interval line w dotted black

    saveas(fig9a,strcat(savePath,'WT_RevLog.png')) % Save displayed fig img
%     saveas(fig9a,strcat(savePath,'WT_RevTr9.png')) % Use for linear figur

%%
% Run this section after running the first section to load your variables

close all % Close any open figures
clc % Clear command window

% Repeat this code for as many .mat variables you want to load
loaded = load('Y:\Lena\Analyses\ExtraAnalyses\.mat Variables\WT_edges.mat');
variableName = loaded.variableNameInMatFile; % Redefine name of variable
isNZ_variableName{1} = (~variableName{1}==0); % Replaces 0's w/ NaN, left C
isNZ_variableName{2} = (~variableName{2}==0); % Replaces 0's w/ NaN, right

% Repeat below code for as many .mat variables you want statistics for
big = vertcat(variableName{1}(isNZ_variableName{1}),variableName{2}(isNZ_variableName{2}));
variMean = mean(big); % Gets mean (average) of left and right "cleaned" val
varistdev = std(big); % Gets standard deviation of left/right "cleaned" val
varin = length(big); % Gets # of cells represented in data (left+right Cs)

position = 1;
rowLeft = variableName{1}(position,:); % Getting all left Tr9 cell data
rowRight = variableName{2}(position,:); % Getting all right Tr9 cell data
big = horzcat(rownode1(isNZ_nodePlotPlot{1}(position,:)),rownode2(isNZ_nodePlotPlot{2}(position,:)));
variTr9Mean = mean(big); % Gets mean (average) of left and right "cleaned"
variTr9stdev = std(big); % Gets standard deviation of left/right "cleaned"
variTr9n = length(big); % Gets # of cells represented in data (left+right)
