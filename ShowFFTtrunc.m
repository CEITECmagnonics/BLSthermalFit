
Fskx = 2*pi/StepX;
Fsky = 2*pi/StepY;

kx = (-sizeOfData(1)/2:sizeOfData(1)/2-1)*(Fskx)/sizeOfData(1);
ky = (-sizeOfData(2)/2:sizeOfData(2)/2-1)*(Fsky)/sizeOfData(2);

figure('name', 'Ky=0');
[XF, FX] = ndgrid(kx,fS/1e9);
hKspace = surf(XF,FX,abs(squeeze(FFTtruncConj(:,end/2,:))));
hold on
plot3(kPython/1e6, fBVpython,linspace(1e40,1e40,501));
plot3(kPython/1e6, fBV11python,linspace(1e40,1e40,501));
set(hKspace,'edgecolor','none')
xlabel('Kx (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 50, 4 20]);
title('Dispersion Ky=0');
set(gca,'FontSize',30)

figure('name', 'Density of states - python')
plot(fS/1e9, sum(squeeze(FFTtruncConj(:,end/2,:)),1)./max(max(sum(squeeze(FFTtruncConj(:,end/2,:)),1))), fBVpython(1:end-1), abs(fBVDoSpython)./max(abs(fBVDoSpython)))
xlabel('f (GHz)')
ylabel('Density of states')

figure('name', 'Kx=0');
[YF, FY] = ndgrid(ky,fS/1e9);
hKspace = surf(YF,FY,(abs(squeeze(FFTtruncConj(end/2,:,:)))));
hold on;
plot3(kPython/1e6, fDEPython,linspace(1e40,1e40,501));
plot3(kPython/1e6, FDE11Python,linspace(1e40,1e40,501));
set(hKspace,'edgecolor','none')
xlabel('Ky (rad/\mum)');
ylabel('f (GHz)');
view(2);
axis([0 50, 0 20]);
title('Dispersion Kx=0');
set(gca,'FontSize',30)

figure('name', 'Density of states - python')
plot(fS/1e9, sum(squeeze(FFTtruncConj(end/2,:,:)),1)./max(max(sum(squeeze(FFTtruncConj(end/2,:,:)),1))), fDEPython(1:end-1), fDEDoSpython./max(fDEDoSpython))
xlabel('f (GHz)')
ylabel('Density of states')