% This function creates a 2 dimentional Hann window

%   Example: 
%   To compute windowed 2D fFT
%   [r,c]=size(img);
%   w=window2(r,c,@hamming);
% 	fft2(img.*w);

function w=hann(N,M)
for i = 1:N
    wc(i)=0.5*(1-cos(2*pi*i/N));
end
for i = 1:M
    wr(i)=0.5*(1-cos(2*pi*i/M));
end
[maskr,maskc]=meshgrid(wr,wc);

%maskc=repmat(wc,1,M); Old version
%maskr=repmat(wr',N,1);

w=maskr.*maskc;

end

