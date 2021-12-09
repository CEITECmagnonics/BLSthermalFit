function [DoS] = fitBLSthermalLumerical(f,A,dc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global fS FFTtrunc BLSkDet
% 
% [KX, KY] = meshgrid(kx, ky);
% dimX = sizeOfData(1);
% dimY = sizeOfData(2);
% FsCellX = (2*pi)/StepX;
% FsCellY = (2*pi)/StepY;
% kxi = (-dimX/2:dimX/2-1)*(FsCellX)/dimX;
% kyi = (-dimY/2:dimY/2-1)*(FsCellY)/dimY;
% [KXi, KYi] = meshgrid(kxi, kyi);
% BLSkDet = interp2(KX, KY, I2, KXi, KYi);

% figure('name', 'Lumerical result')
% surf(KXi, KYi, BLSkDet, 'edgecolor','none')
% view(2)
% xlim([0 120]);
% ylim([0 120]);
% caxis([1e4 1e7]);


% Bose-Einstein distribution
cL = 299792458;
wi = 2*pi*cL/532e-9; % Frequency of light
hbar = 1.0545718e-34;
kb = 1.38064852e-23;
T = 300;
BE = wi*(wi+2*pi*abs(fS)).^3.*(1./(exp(hbar*2*pi*abs(fS)/(kb*T))-1))./1e69;
BE = (BE./max(abs(BE)))';
% plot(abs(fS), BE)

% DoS
DoSsim = squeeze(sum(sum(abs(BLSkDet.*(FFTtrunc.*conj(FFTtrunc))/max(max(max(FFTtrunc.*conj(FFTtrunc))))),2),1));
DoSsim = DoSsim.*BE./max(DoSsim.*BE);
DoS = double(A*interpn(abs(fS)/1e9,DoSsim,f,'linear',0)+dc);
end