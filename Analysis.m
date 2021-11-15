step = [3.75*4e-9, 3.75*4e-9, 3.75*4e-9];

% Data = DataInterp(:,:,:);

sizeOfData = size(Data);
dt = 8e-12; % Time step in seconds
Fs = 1/dt;
StepX = step(1)*1e6;
StepY = step(1)*1e6;
Fskx = 2*pi/StepX;
Fsky = 2*pi/StepY;

figure('name', 'Real space before filtering');
[X,Y] = ndgrid(StepX:StepX:sizeOfData(1)*StepX,StepY:StepY:sizeOfData(2)*StepY);
surf(X,Y,squeeze(abs(Data(:,:,200))));
view(2);
axis 'tight';
axis square;
xlabel('x (\mum)');
ylabel('y (\mum)');
shading interp;
% % ~~~~Xt map
xVsTmean = squeeze(mean(Data(:,:,:),2));
figure('name', 'x vs time before filtering');
x = StepX:StepX:sizeOfData(1)*StepX;
t = dt:dt*1e9:sizeOfData(3)*dt*1e9;
[X,T] = ndgrid(x,t);
surf(X,T,xVsTmean);
view(2);
axis 'tight';
xlabel('x (\mum)');
ylabel('t (ns)');
shading interp;

f = (-sizeOfData(3)/2:sizeOfData(3)/2-1)*(Fs)/sizeOfData(3);
kx = (-sizeOfData(1)/2:sizeOfData(1)/2-1)*(Fskx)/sizeOfData(1);
ky = (-sizeOfData(2)/2:sizeOfData(2)/2-1)*(Fsky)/sizeOfData(2);
[~,kyIndex] = min(abs(ky-0));

fprintf('~~~~~~~ Clean data before FFT ~~~~~~~\n');
tic
[r, c] = size(Data(:,:,1));
win = hann(r,c);
for i = 1:sizeOfData(3)
    Data(:,:,i) = detrend_2d(squeeze(Data(:,:,i))).*win;
end
toc


fprintf('~~~~~~~ Perform FFT ~~~~~~~\n');
tic
FFT = fftshift(fftn(Data));
toc


[~,FreqIndex] = min(abs(f-28e9));
figure('name', 'K-space for one frequency');
[KX, KY] = ndgrid(kx,ky);
hKspace = surf(KX,KY,abs(FFT(:,:,FreqIndex)).^2);
set(hKspace,'edgecolor','none')
xlabel('Kx (rad/\mum)');
ylabel('Ky (rad/\mum)');
view(2);
axis([0 100, 0 100]);
title('f=28 GHz');
set(gca,'FontSize',30)


figure('name', 'Ky=0');
[XF, FX] = ndgrid(kx,f/1e9);
hKspace = surf(XF,FX,(abs(squeeze(FFT(:,end/2,:)))).^2);
hold on
plot3(kPython, fBVPython,linspace(1e40,1e40,150));
plot3(kPython, f11,linspace(1e40,1e40,150));
set(hKspace,'edgecolor','none')
xlabel('Kx (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 100, 20 40]);
title('Dispersion Ky=0');
set(gca,'FontSize',30)

figure('name', 'Kx=0');
[YF, FY] = ndgrid(ky,f/1e9);
hKspace = surf(YF,FY,(abs(squeeze(FFT(end/2,:,:)))).^2);
hold on;
plot3(kPython, fDePython,linspace(1e40,1e40,150));
plot3(kPython, f11DE,linspace(1e40,1e40,150));
set(hKspace,'edgecolor','none')
xlabel('Ky (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 100, 20 40]);
title('Dispersion Kx=0');
set(gca,'FontSize',30)

% FMR
mu0 = 4*pi*1e-7;
Ms = 800e3; %A/m
Hext = 550*1e-3/mu0; %A/m
d = 30e-9; %m
gamma = 28.8; %GHz/T
Aex = 16e-12; %J/m

w0 = mu0*gamma*Hext;
wM = mu0*gamma*Ms;
A = Aex*2/(Ms^2*mu0); % 
weff = 3.0e-6;
FMR = sqrt(w0*(w0 + wM));

%BLS k detection
BLSkDet = BLSkDetectionGauss(StepX, StepY, sizeOfData(1), sizeOfData(2), 4.3, 1);
figure('name', 'Detection function');
[KX, KY] = ndgrid(kx,ky);
hKspace = surf(KX,KY,BLSkDet);
set(hKspace,'edgecolor','none')
xlabel('Kx (rad/\mum)');
ylabel('Ky (rad/\mum)');
view(2);
axis([0 30, 0 30]);
title('Detection function');
set(gca,'FontSize',30)


figure('name', 'Kx=0 BLS detection');
[YF, FY] = ndgrid(ky,f/1e9);
BLSconvFFTky = FFT.*BLSkDet;
hKspace = surf(YF,FY,(abs(squeeze(BLSconvFFTky(end/2,:,:)))).^2);
hold on;
plot3(kPython, fDePython,linspace(1e40,1e40,150));
plot3(kPython, f11DE,linspace(1e40,1e40,150));
set(hKspace,'edgecolor','none')
xlabel('Ky (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 50, 20 40]);
title('Dispersion Kx=0');
set(gca,'FontSize',30)


% Trunctated simulated dispersion to only fundamental mode
[~,FreqIndexMax] = min(abs(abs(f(1:round(end/2)))-60e9));
[~,FreqIndexMin] = min(abs(abs(f(1:round(end/2)))-20e9));
FFTtrunc = FFT(:,:,FreqIndexMax:FreqIndexMin);
fS = -f(FreqIndexMax:FreqIndexMin);
FFTtrunc(:,:,1:FreqIndexMax) = 0;
% FFTtrunc(:,:,end-FreqIndexMax:end) = 0;



% Density of states
figure('name', 'DoS');
DoS = squeeze(sum(sum(abs(BLSkDet.*FFT.^2),2),1));
cL = 299792458;
wi = 2*pi*cL/532e-9; % Frequency of light
hbar = 1.0545718e-34;
kb = 1.38064852e-23;
T = 300;
BE = wi*(wi+2*pi*f).^3.*(1./(exp(hbar*2*pi*f/(kb*T))-1))./1e46;
BE = BE';
plot(DoS./max(DoS),f/1e9, DoS.*BE./max(DoS.*BE), f/1e9, [0 1], [FMR FMR], 'b--', (BLSexp)./max(BLSexp), fexp, '-o')
xlabel('Density of states ()');
ylabel('Frequency (GHz)');
legend('DoS', 'DoS with Bose-Einstein', 'FMR', 'Experimental data')
ylim([20 45]);
xlim([0 1]);
global sizeOfData StepX StepY fS FFTtrunc