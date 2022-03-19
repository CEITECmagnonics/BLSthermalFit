function [BLSMatrix] = BLSkDetectionGauss(fx, fy, sigX,sigY, A)
%Compute matrix of BLS k detection, !steps in micrometer!
%Prepare matrix for cyclus
dimX = length(fx);
dimY = length(fy);
kStepX = fx(2) - fx(1);
kStepY = fy(2) - fy(1);
BLSMatrix = zeros(length(fx),length(fy));
%Compute center of matrix
sx = round((dimX+1)/2); sy = round((dimY+1)/2);
% sig = 4.3; %Detection limit of BLS in rad/um
% A = 1;
for i=1:length(fx)
    for j=1:length(fy)
            BLSMatrix(i,j) = A*exp(-((((i-sx)*kStepX).^2)/(2*sigY.^2)))*exp(-((((j-sy)*kStepY).^2)/(2*sigX.^2))); 
    end
end