%Homework 4 part 2
%By Yejing Zhou & Mengyan Li

%plant model
% Ap = [-1.3046, 1.0, -0.21420, 0; 47.411, 0, -104.83, 0; 0, 0, 0, 1.0; 0, 0, -1.2769, -1.356];
% Bp = [0;0;0;12769.0];

clear
clc
% actuator model
wn = 2*pi*11;
xi = 0.707;
% plant model with actuator
Ap = [-1.3046,1.0, -0.21420, 0; 47.711, 0, -104.83, 0;0, 0, 0, 1; 0, 0, -wn^2, -2*xi*wn];
Bp = [0;0;0;wn^2];
Cp = [-1156.9, 0, -189.95, 0; eye(4)];
Dp = 0.*Cp*Bp;
%plant model data
Za = -1156.9;
Zs = -189.95;
V = Za/(-1.3046);
Ma = 47.711;
Ms = -104.83;
% RSLQR design zdot = Arz + Brn
%state feedback design
Ar = [0, Za, 0; 0, Za/V, 1; 0, Ma, 0];
Br = [Zs;Zs/V;Ms];

% Solve ARE
Q = 0.*Ar;
R = 1;
xeig = [];
qq = logspace(-6,-1,50);
xol = eig(Ar); 

% set scales
t = 0:0.002:2;
w = logspace(-1,3,500);
% allocate matrix
Az_st = 0.*ones(numel(qq), numel(t));%acceleration
q_st = 0.*ones(numel(qq), numel(t));%pitch rate
Del_st = 0.*ones(numel(qq), numel(t)); %elevon
Deldot_st = 0.*ones(numel(qq), numel(t)); %differentila elevon


% calculate gain values with different q11 values.
    Q(1,1) = qq(20);
    [Kx_lqr, ~, ~] = lqr(Ar, Br, Q, R);
   
    % populate the controller matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0. ];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;
    
    % connect the controller with plant, form a closed loop system
    z = inv(eye(size(Dc1*Dp)) - Dc1*Dp);
    Acl = [(Ap + Bp*z*Dc1*Cp) (Bp*z*Cc); (Bc1*(Cp+Dp*z*Dc1*Cp)) (Ac+Bc1*Dp*z*Cc)];
    Bcl = [Bp*z*Dc2; (Bc2+Bc1*Dp*z*Dc2)];
    Ccl = [(Cp+Dp*z*Dc1*Cp) (Dp*z*Cc)];
    Dcl = (Dp*z*Dc2);
    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
   
    %compute the eigenvalues for Acl to plot root locus plot
    xacl = eig(Acl);
    
    %loop gain model at input breakpoint
    Ali = [Ap, 0.*Bp*Cc; Bc1*Cp, Ac];
    Bli = [Bp; Bc1*Dp];
    Cli = -[Dc1*Cp, Cc];
    Dli = -[Dc1*Dp];
    sys_li = ss(Ali,Bli,Cli,Dli);
    
    %loop gain model at output breakpoint
    Alo = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
    Blo = [ Bp*Dc1; Bc1];
    Clo = -[ Cp Dp*Cc];%change sign for loop gain
    Dlo = -[ Dp*Dc1];
    sys_lo = ss(Alo,Blo,Clo,Dlo);
    
    %frequency analysis
    Li = freqresp(sys_li,w);
    Lo = freqresp(sys_lo,w);
    T = freqresp(sys_cl,w);
    S = 1 - T;
    
    sr = sigma(sys_li,w,3);
    sr_min = min(abs(sr)); 
    rd = sigma(sys_li,w,2);
    rd_min = min(abs(rd)); 
   
    
    [nCp,nAp] = size(Cp);
    [~,nBp]   = size(Bp);
    
    theta = -2*pi:.001:2*pi;
    x1 = cos(theta);
    y1 = sin(theta);
    %analysis at input breakpoint
    for i = 1:numel(w),
    s = sqrt(-1)*w(i);
    G = Cp*inv(s*eye(size(Ap)) - Ap)*Bp + Dp;
    K = Cc*inv(s*eye(size(Ac)) - Ac)*Bc1 + Dc1;
    Lof(i) = -K*G;
    rd(i) = 1.+ Lof(i);
    sr(i) = 1.+1./Lof(i);
    sysin = Cli*inv(s*eye(size(Ali)) - Ali)*Bli + Dli;
    sysout = Clo*inv(s*eye(size(Alo)) - Alo)*Blo + Dlo;
    
    for j = 1:nBp,
        F = eye(size(sysin));
        F(j,j) = 0.;
        Tin(j,i) = inv(eye(size(sysin)) + sysin*F)*sysin;
    end
    
    for j = 1:nCp,
        F = eye(size(sysout));
        F(j,j) = 0.;
        Tout0 = inv(eye(size(sysout)) + sysout*F)*sysout;
        Tout(j,i) = Tout0(j,j);
    end
end
   
% bode plot at input breakpoint
figure
semilogx(w, 20*log10(abs(squeeze(L))),'LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')


%nyquist at input
xu = x1-1;
yu = y1;
xi = rd_min*x1 - 1;
yi = rd_min*y1 ;
figure
plot(xi,yi,'r',xu,yu,'k',real(squeeze(Li)),imag(squeeze(Li)),'b');
grid
axis([-3 2 -2 2]);
legend('Loop','Unit Circle','Li');
xlabel('Real')
ylabel('Imaginary')
title('Nyquist Plot at Plant Input Breakpoint')


% return difference 1+L
figure
semilogx(w,20*log10(abs(rd)),'r');
grid
xlabel('Frequency')
ylabel('Mag dB')
title('Return Difference at Plant Input Breakpoint')

% stability robustness 1+L-1
figure
semilogx(w,20*log10(abs(sr)),'b');
grid
xlabel('Frequency')
ylabel('Mag dB')
title('Stability Robustness at Plant Input Breakpoint')


% calculate the phase margin and gain margin with minimum singular values.
disp('******SV Margins******')
rd_ngm = 1/(1+rd_min);
rd_pgm = 1/(1-rd_min);
rd_pha = 2*asin(rd_min/2);
rd_ngm_dB = 20*log10(rd_ngm)
rd_pgm_dB = 20*log10(rd_pgm)
rd_pha_deg = 180*rd_pha

sr_ngm= 1-sr_min;
sr_pgm= 1+sr_min;
sr_pha = 2*asin(sr_min/2);
sr_ngm_dB = 20*log10(sr_ngm)
sr_pgm_dB = 20*log10(sr_ngm)
sr_pha_deg = 180*sr_pha/pi 

% sensitivity E(s)/R(s) = Scl(s)
figure
for j = 1:2,
    subplot(2,1,j)
    semilogx(w,20*log10(abs(squeeze(S(j,1,:)))),'b');
    grid
    xlabel('Frequency ')
    ylabel('Mag (dB)')
    title('Sensitivity S(s)')
end

% Complememtary sensitivity Y(s)/R(s) = Lcl(s)
figure
for j = 1:2,
    subplot(2,1,j)
    semilogx(w,20*log10(abs(squeeze(T(j,1,:)))),'b');
    grid
    xlabel('Frequency ')
    ylabel('Mag (dB)')
    title('Complementary Sensitivity T(s)')
end

   
