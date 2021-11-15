function [BLSMatrix] = BLSkDetectionGauss(stepX,stepY,dimX,dimY, sig, A)
%Compute matrix of BLS k detection, !steps in micrometer!
%Compute length of steps
FsCellX = (2*pi)/stepX;
FsCellY = (2*pi)/stepY;
fx = (-dimX/2:dimX/2-1)*(FsCellX)/dimX;
kStepX = fx(2)-fx(1);
fy = (-dimY/2:dimY/2-1)*(FsCellY)/dimY;
kStepY = fy(2)-fy(1);
%Prepare matrix for cyclus
BLSMatrix = zeros(length(fx),length(fy));
%Compute center of matrix
sx = round((dimX+1)/2); sy = round((dimY+1)/2);
% sig = 4.3; %Detection limit of BLS in rad/um
% A = 1;
for i=1:length(fx)
    for j=1:length(fy)
            BLSMatrix(i,j) = A*exp(-((((i-sx)*kStepX).^2)/(2*sig.^2)))*exp(-((((j-sy)*kStepY).^2)/(2*sig.^2))); 
    end
end