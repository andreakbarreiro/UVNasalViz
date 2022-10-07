function [prefPhase,valAtPrefPhase] = calcPrefPhase(avgW,timeMod)
% calcPrefPhase: calculate preferred phase at which avgW is greatest for
%   every point in space
%   INPUTS:
%       avgW:      (nS x nT)  quantity
%       timeMod    (nT x 1)   time (or phase) values
%   OUTPUTS:
%       
% Maximum over time

maxOT       = max(avgW,[],2);
prefTimeInd = zeros(size(maxOT));
[nFace,nT]  = size(avgW);

for j1=1:nFace
   prefTimeInd(j1) = find(avgW(j1,:)==maxOT(j1),1);
end

prefPhase   = timeMod(prefTimeInd);

valAtPrefPhase = maxOT;

end

