function [Results] = ChannelCouplingAnalysis(coupling_type, method, deadtime)

%Copyright (c) 2016 by Gary Iacobucci
%
%%%%%%%INPUT VARIABLES:
%coupling_type - specify whether you are analyzing channels that are
%positively or negatively coupled. Enter 1 for positive coupling and 2 for
%negative coupling

%method - specify the method of analyzing the transition probability matrix
%(eg. by sequence of state transitions (where all diagonal elements = 0) or
%by time-series (state at each sampled time point)). Enter 1 for time-series
%method or enter 2 for state-sequence method.

%deadtime - enter the deadtime in the form of sample point number. must be
%a whole integer. a value of 1 implies no imposed deadtime. 

%%%%%%  OUTPUT VARIABLES:
%Results - will contain the results of the fitted parameters from the
%'channelcoupling.m' script for each file. Fitted parameters are based on
%the notation from Chung and Kennedy, 1996 and Moreno et al. 2016 eLife
% output.k - coupling coefficient
% output.r - open-to-open probability
% output.z - closed-to-closed probability


%%%%%  HOW TO USE
%prior to running this script, you will need to create an array named 'fileID' containing
%the file names of the files you will analyze. This program is written to
%accept files in the DWT format from QUB software that expresses state
%index and corresponding dwell time in column format. 
%
%example array format:
%
% fileID = {'GI-20160810a.dwt';
% 'file1.dwt';
% 'file2.dwt';
% 'file3.dwt';
% 'file4.dwt';
% 'file5.dwt'};
%
%Example command line input to run script:
% [Results] = ChannelCouplingAnalysis(2, 1, 1);


Results = {};
h = waitbar(0,'Please wait...');
for n = 1:1:numel(fileID)
    waitbar(n/numel(fileID))
    filename = fileID{n,1};
    [output]=channelcoupling(filename,coupling_type,method,deadtime);
    Results{n,1} = output;
end
close(h)