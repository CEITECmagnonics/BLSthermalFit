function [DoS] = fitBLSthermal(f,FWTMx, FWTMy, A,dc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global sizeOfData StepX StepY fS FFTtrunc

sigX = FWTMx/4.29193;
sigY = FWTMy/4.29193;
BLSkDet = BLSkDetectionGauss(StepX, StepY, sizeOfData(1), sizeOfData(2), sigX, sigY, 1);

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