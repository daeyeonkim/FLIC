% Analysis of data from FLIC
% Developed by Daeyeon Kim
% Created on Mar 09, 2020
% Last updated date: April 23, 2020

%% clean up
close all
clear

%% set variables
analTime = 1440; % in min
tStep = 0.2; 
maxFrame = 60*analTime/tStep; % maximum frame to be analyzed
%onSetFrame = 10*60/tStep;
onSetFrame = 1;
windowSize = 5*60/tStep; % window size for calculating baseline (frame)
thRemoveChannel = 50; % threshold for removing channel data in baseline substracted signal value
thContactPeak = 10; % threshold for removing irrelanvent signal for contact
% index of the channel to remove
idxRemoveChannel = [];
%idxRemoveChannel = find(meanBaseLine > thRemoveChannel);
GROUP = {'CS'};
flagTwoChoiceAssay = 1; % option for calculating preference idex (0: OFF, 1: ON)
exportFlag = 1; % option for exporting the result (0: OFF, 1: ON)

%% read input file
[filename,path] = uigetfile('*.csv');
if isequal(filename,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,filename)]);
end

% import the data as a table
flicTable =  readtable(strcat(path,filename),...
    'Delimiter',',','ReadVariableNames',true,'HeaderLines',0,'Format','%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f');

% check the maximum frame
if size(flicTable,1) < maxFrame
    maxFrame = size(flicTable,1);
end

% get time vector
timeVec = 0.2/60*(onSetFrame:maxFrame); % in min

% rearrange the data as a table
flicTable = flicTable(onSetFrame:maxFrame,:);
flicTable = addvars(flicTable,timeVec','Before','W1','NewVariableNames','ElapsedTime');

% get variables in the table
vars = flicTable.Properties.VariableNames(6:end);

%% analize the flic signal
% get the flic data
flicData = table2array(flicTable(:,6:end));
% calculate the baseline
baseLine = movmedian(flicData,windowSize,1);
% calculate the mean baseline
meanBaseLine = mean(baseLine,1);

% plot the flicData to identify channels to remove in the anlaysis
figure(1)
s = stackedplot(flicTable,vars,'XVariable','ElapsedTime','Marker','.');
myColours = [0,0,1;1,0,0]; % this is your 12x3 matrix of colours (one for each line)
for k = 1:length(s.LineProperties)
        switch mod(k,2)
            case 1
                s.LineProperties(k).Color = myColours(1,:);
                s.LineProperties(k).MarkerFaceColor = myColours(1,:);
                s.LineProperties(k).MarkerEdgeColor = myColours(1,:);
                s.AxesProperties(k).YLimits = [0 500];
                %s.LineProperties(k).PlotType = 'scatter';
            case 0
                s.LineProperties(k).Color = myColours(2,:);
                s.LineProperties(k).MarkerFaceColor = myColours(2,:);
                s.LineProperties(k).MarkerEdgeColor = myColours(2,:);
                s.AxesProperties(k).YLimits = [0 500];
        end        
end
set(gca,'FontSize',18)

% remove the unfavorable channel data 
% index of the channel to be analyzed
analChannel = setdiff(1:numel(vars),idxRemoveChannel);
% calculate the filc signal
flicSignal = flicData - baseLine;
% thresholding the flicSignal
flicLabel = flicSignal > thContactPeak;

% set variables 
durContact = []; % duration of contact
peakSignal = []; % peak value of contact
meanDurContact = []; % mean duraion of contact
meanIntervalContact = []; % mean interval of contacts
totDurContact = []; % total duration of all contacts 
totDurContactChannel = nan(1,12); % total duration of contact in channel
channelA = 1:2:11;
channelB = 2:2:12;
groupIdx = [];

for i = 1:numel(analChannel)
    tempLabel = flicLabel(:,analChannel(i));
    [start,stop] = segmentBehavior(tempLabel',maxFrame);
    tempDurContact = stop - start;
    tempIntervalContact = start(2:end) - stop(1:end-1);
    tempIntervalContactMin = tempIntervalContact*tStep/60;
    tempDurContactSec = tempDurContact*tStep
    tempTotDurContactFly = nansum(tempDurContactSec);
    tempMeanDurContactFly = nanmean(tempDurContactSec);
    tempMeanIntervalContactFly = nanmean(tempIntervalContactMin);
    durContact = cat(2,durContact,tempDurContactSec);
    totDurContact = cat(2,totDurContact,tempTotDurContactFly);
    meanDurContact = cat(2,meanDurContact,tempMeanDurContactFly);
    meanIntervalContact = cat(2,meanIntervalContact,tempMeanIntervalContactFly);
      
    tempPeak = [];
    tempSignal = flicSignal(:,analChannel(i));
    for j=1:numel(start)
        tempPeak(j) = max(tempSignal(start(j):stop(j)));
    end
    peakSignal = cat(2,peakSignal,tempPeak);
end

% assign group for the box plot
groupIdx = cat(1,groupIdx,ones(numel(meanDurContact),1));
sigma = 0.05;

% plot results
figure(2)
subplot(2,4,1)
boxplot(meanDurContact,groupIdx,'Labels',GROUP)
hold on
scatter(1.25+sigma*randn(numel(meanDurContact),1),meanDurContact);
set(gca,'TickDir','out','TickLength',[0.005 0.005]); % The only other option is 'in'
set(gca,'FontSize',12);
ylabel('Duration of contact (sec)')

subplot(2,4,2)
boxplot(meanIntervalContact,groupIdx,'Labels',GROUP)
hold on
scatter(1.25+sigma*randn(numel(meanIntervalContact),1),meanIntervalContact);
set(gca,'TickDir','out','TickLength',[0.005 0.005]); % The only other option is 'in'
set(gca,'FontSize',12);
ylabel('Interval of contact (min)')

if flagTwoChoiceAssay
    totDurContactChannel(analChannel) = totDurContact;
    totDurContactChannelA = totDurContactChannel(channelA);
    totDurContactChannelB = totDurContactChannel(channelB);
    preferIdx = (totDurContactChannelA - totDurContactChannelB)./(totDurContactChannelA+totDurContactChannelB);
    
    subplot(2,4,3)
    boxplot(preferIdx,'Labels',GROUP)
    hold on
    scatter(1.25+sigma*randn(numel(preferIdx),1),preferIdx);
    set(gca,'TickDir','out','TickLength',[0.005 0.005]); % The only other option is 'in'
    set(gca,'FontSize',12);
    ylabel('Preference index')
    
end


if exportFlag == 1
% --------------------------------------------------------
% export the results
% switching frequency in time section
fileNameExcel = strcat(path,'flicData_',filename,'.xlsx');
col_header={'Duration of contact (sec)','Interval of contact (min)'};     %Row cell array (for column labels)
row_header={'Channel'};     %Column cell array (for row labels)

dataAll = [cat(1,meanDurContact',nanmedian(meanDurContact)), cat(1,meanIntervalContact',nanmedian(meanIntervalContact))];

xlswrite(fileNameExcel,dataAll,'Sheet1','B2');     %Write data
xlswrite(fileNameExcel,col_header,'Sheet1','B1');     %Write column header
xlswrite(fileNameExcel,row_header,'Sheet1','A1');      %Write row header
xlswrite(fileNameExcel,cat(1,num2cell(analChannel'),'median'),'Sheet1','A2');

col_header={'Preference Index','total duration of contact in the channel A (sec)', 'total duration of contact in the channel B (sec)'};     %Row cell array (for column labels)
row_header={'Chamber'};     %Column cell array (for row labels)

chamberIdx = find(~isnan(preferIdx));
dataPI = [cat(1,preferIdx(chamberIdx)',nanmedian(preferIdx(chamberIdx))), cat(1,totDurContactChannelA(chamberIdx)',nanmedian(totDurContactChannelA(chamberIdx))),...
    cat(1,totDurContactChannelB(chamberIdx)',nanmedian(totDurContactChannelB(chamberIdx)))];

xlswrite(fileNameExcel,dataPI,'Sheet2','B2');     %Write data
xlswrite(fileNameExcel,col_header,'Sheet2','B1');     %Write column header
xlswrite(fileNameExcel,row_header,'Sheet2','A1');      %Write row header
xlswrite(fileNameExcel,cat(1,num2cell(chamberIdx'),'median'),'Sheet2','A2');
end

meanPI = nanmean(preferIdx)
medianPI = nanmedian(preferIdx)
preferIdx

% calculate the CDF
binSize = 0.6;
edges = 0:binSize:200;
histDurContact = histcounts(durContact,edges);
pdfDurContact = histDurContact/sum(histDurContact)/binSize;
estimPdfDurContact = ksdensity(durContact,edges,'Bandwidth',1);
%cdfDurContact = cumsum(pdfDurContact)*binSize;

function [start,stop] = segmentBehavior(behavior,maxFrame)
partition_idx = [];
partition_idx = [0 behavior 0];
start = find(~partition_idx(1:(end-1)) &  partition_idx(2:end));
stop = find( partition_idx(1:(end-1)) & ~partition_idx(2:end));

% truckate frames out of the total analyzing time
trunkIdx = find(stop > maxFrame,1);
if ~isempty(trunkIdx)
    switch start(trunkIdx) > maxFrame
        case 0
            stop(trunkIdx) = maxFrame;
            start(trunkIdx+1:end) = [];
            stop(trunkIdx+1:end) = [];
        case 1
            start(trunkIdx:end) = [];
            stop(trunkIdx:end) = [];
    end
end
end