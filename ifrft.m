function ifaf = ifrft(Faf, a)
% The fast Inverse Fractional Fourier Transform
% input: Faf = Fractional Fourier coefficients
%        a = fractional power
% output: ifaf = fast Inverse Fractional Fourier transform

% H.M. Ozaktas, M.A. Kutay, and G. Bozdagi.
% [i]Digital computation of the fractional Fourier transform.[/i]
% IEEE Trans. Sig. Proc., 44:2141--2150, 1996.

error(nargchk(2, 2, nargin));
Faf = Faf(:);
N = length(Faf);
shft = rem((0:N-1)+fix(N/2),N)+1;
sN = sqrt(N);
a = mod(a,4);
% do special cases
if (a==0), ifaf = Faf; return; end
if (a==2), ifaf = flipud(Faf); return; end
if (a==1), ifaf(shft,1) = ifft(Faf(shft))*sN; return; end 
if (a==3), ifaf(shft,1) = fft(Faf(shft))/sN; return; end
% reduce to interval 0.5 < a < 1.5
if (a>2.0), a = a-2; Faf = flipud(Faf); end
if (a>1.5), a = a-1; Faf(shft,1) = ifft(Faf(shft))*sN; end
if (a<0.5), a = a+1; Faf(shft,1) = fft(Faf(shft))/sN; end
% the general case for 0.5 < a < 1.5
alpha = a*pi/2;
tana2 = tan(alpha/2);
sina = sin(alpha);
% chirp post multiplication
ifaf = exp(1i*(1-a)*pi/4)*Faf(N:2:end-N+1);
% inverse chirp convolution
c = pi/N/sina/4;
ifaf = fconv(exp(-1i*c*(-(4*N-4):4*N-4)'.^2),ifaf);
%%%%%%%%%%   ifaf = ifaf(4*N-3:8*N-7)*sqrt(c/pi);
ifaf = ifaf(2*N-1:5*N-3)*sqrt(c/pi);
% inverse chirp premultiplication
chrp = exp(1i*pi/N*tana2/4*(-2*N+2:2*N-2)'.^2);
%%%%%%%%  ifaf = chrp.*ifaf
ifaf = chrp(1:N).*ifaf(1:N);

% sinc interpolation
ifaf = interp(ifaf);
end

function xint=interp(x)
% sinc interpolation
N = length(x);
y = zeros(2*N-1,1);
y(1:2:2*N-1) = x;
xint = fconv(y(1:2*N-1), sinc((-(2*N-3):(2*N-3))'/2));
xint = xint(2*N-2:end-2*N+3);
end

function z = fconv(x,y)
% convolution by fft
N = length([x(:);y(:)])-1;
P = 2^nextpow2(N);
z = ifft( fft(x,P) .* fft(y,P));
z = z(1:N);
end
