
Fskx = 2*pi/StepX;
Fsky = 2*pi/StepY;

kx = (-sizeOfData(1)/2:sizeOfData(1)/2-1)*(Fskx)/sizeOfData(1);
ky = (-sizeOfData(2)/2:sizeOfData(2)/2-1)*(Fsky)/sizeOfData(2);

figure('name', 'Ky=0');
[XF, FX] = ndgrid(kx,fS/1e9);
hKspace = surf(XF,FX,log(abs(squeeze(FFTtrunc(:,end/2,:)))).^2);
hold on
plot3(kPython, fBVPython,linspace(1e40,1e40,150));
plot3(kPython, f11,linspace(1e40,1e40,150));
set(hKspace,'edgecolor','none')
xlabel('Kx (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 150, 20 50]);
title('Dispersion Ky=0');
set(gca,'FontSize',30)

figure('name', 'Kx=0');
[YF, FY] = ndgrid(ky,fS/1e9);
hKspace = surf(YF,FY,log((abs(squeeze(FFTtrunc(end/2,:,:)))).^2));
hold on;
plot3(kPython, fDePython,linspace(1e40,1e40,150));
plot3(kPython, f11DE,linspace(1e40,1e40,150));
set(hKspace,'edgecolor','none')
xlabel('Ky (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 150, 20 50]);
title('Dispersion Kx=0');
set(gca,'FontSize',30)