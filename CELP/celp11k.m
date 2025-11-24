function [xhat,e,k,theta0,P,b] = celp11k(x,N,L,M,c,cb,Pidx)
%  celp --> ~11000 bps CELP analyzer and synthesizer.
%
%   [xhat,e,k,theta0,P,b] = celp11k(x,N,L,M,c,cb,Pidx)
%
%   Bits used:
%      LPC coefficients: 8 bits each
%      Gain theta0:       9 bits
%      Pitch filter b:    10 bits
%
%   Total ≈ 220 bits / frame → ~11 kbps @ 20 ms frames

Nx = length(x);                         
F  = fix(Nx/N);                         
J  = N/L;                               

% Initialize output
xhat   = zeros(Nx,1);                   
e      = zeros(Nx,1);                   
k      = zeros(J,F);                    
theta0 = zeros(J,F);                    
P      = zeros(J,F);
b      = zeros(J,F);

ebuf  = zeros(Pidx(2),1);               
ebuf2 = ebuf; 
bbuf = 0;  
Zf = []; Zw = []; Zi = [];

for f = 1:F

    n = (f-1)*N+1:f*N;

    % ---- CELP ANALYZER ----
    [kappa,kf,theta0f,Pf,bf,ebuf,Zf,Zw] = celpana( ...
        x(n), L, M, c, cb, Pidx, bbuf, ebuf, Zf, Zw );

    % ======== 11 kbps QUANTIZATION SCHEME ==========
    % LPC coeffs (kappa) – 8 bits each
    sigma = 2/pi*asin(kappa);
    sigma = udecode(uencode(sigma,8),8);
    kappa = sin(pi/2*sigma);

    % Gain θ0 – 9 bits
    theta0 = udecode(uencode(theta0,9,0.2),9,0.2);

    % Pitch filter coefficient b – 10 bits
    b = udecode(uencode(b,10,1.4),10,1.4);
    % ===============================================

    % ---- CELP SYNTHESIS ----
    [xhat(n),ebuf2,Zi] = celpsyn(cb,kappa,kf,theta0f,Pf,bf,ebuf2,Zi);

    % Output excitation & params
    e(n)        = ebuf(Pidx(2)-N+1:Pidx(2));
    k(:,f)      = kf;
    theta0(:,f) = theta0f;
    P(:,f)      = Pf;
    b(:,f)      = bf;

    bbuf = bf(J);       
end
