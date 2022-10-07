function [avgW,timeMod] = calcAvgOverCycles(warray,tarray,Tperiod)
% calcAvgOverCycles : calculate average over cycles (such as a breath
%       cycle)
%    INPUTS:   
%           warray:   (nS x nT) Quantity to average
%           tarray:   (nTx1)    Time steps: ASSUME EQUALLY SPACED!!!
%           Tperiod:  Period
% 
%    OUTPUT
%           avgW:     (nS x nP) Averaged
%           timeMod:  (nPx1)    (Time mod period)
%   

% Check orientation of warray
[m,n]=size(warray);  nT = length(tarray);
if ((m~=nT) && (n~=nT))
    error('warray is not commensurate with tarray');
elseif (n~=nT)  %We want the SECOND dimension to match time
    warray = warray';
end

timeMod     = mod(tarray,Tperiod);

% How many periods do we have?
allFirst = find(abs(timeMod-timeMod(1))<1e-5);
nseq1     = length(allFirst);
onePeriod = timeMod(1:(allFirst(2)-1));

% This is a little baroque
nseq2=length(find(abs(timeMod-onePeriod(end))<1e-5));


% [nseq1 nseq2]

if (nseq1~=nseq2)
    warning('tarray does not contain an integer # of periods: truncating data');
    warray(:,allFirst(end):end) =   [];
end

oldDim = size(warray);
newDim = [oldDim(1),oldDim(2)/nseq2,nseq2];

% oldDim
% newDim(2:end)

warray_reshape = reshape(warray,newDim);

%Average across periods
avgW    = mean(warray_reshape,3);

%Remove extra entries of timeMod
timeMod(allFirst(2):end)=[]; 
end

