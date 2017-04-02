
%*********************
clear
clc
%*******************
%scale
w=logspace(-3,3,1000);
t = 0:0.002:2;
qq = logspace(-6,-1,50);
rtd = 180/pi;
% actuator model
wn = 2*pi*11;
xi = 0.707;
%plant model parameters
Za = -1156.9;
Zs = -189.95;
V = Za/(-1.3046);
Ma = 47.711;
Ms = -104.83;
%plant model 
Ap = [Za/V, Za, 0, Zs;Ma/Za, 0, Ms-(Ma*Zs)/Za, 0;0, 0, 0, 1;0, 0, -wn^2, -2*xi*wn];
Bp = [0; 0; 0; wn^2];
Cp = eye(4);
Dp = 0.*Cp*Bp;
%state feedback model
Ar = [0, 1, 0, 0, 0; 0, Za/V, Za, 0, Zs;0, Ma/Za, 0, Ms-(Ma*Zs)/Za, 0;0, 0, 0, 0, 1;0, 0, 0, -wn^2, -2*xi*wn];
Br = [0; 0; 0; 0; wn^2];

%solve eigenstructure.
Q = 0.*Ar;
R = 1;
xeig = [];
xol = eig(Ar); 
    
Q(1,1) = 1.898e-04;

[Kc,Pr] = lqr(Ar,Br,Q,R);


% populate the controller matrices
Ac =  0.;
Bc1 = [1. 0. 0. 0.];
Bc2 =  -1;
Cc = -Kc(1);
Dc1 = -Kc(2:5);
Dc2 = 0.;

% 5-state RSLQR Frequency Analysis
Ain_lqr = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin_lqr = [ Bp; Bc1*Dp];
Cin_lqr = -[ Dc1*Cp Cc];%change sign for loop gain
Din_lqr = -[ Dc1*Dp];

Lin_lqr = ss(Ain_lqr,Bin_lqr,Cin_lqr,Din_lqr);

%SS model of loop gain at the plant output
Aout_lqr = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
Bout_lqr = [ Bp*Dc1; Bc1];
Cout_lqr = -[ Cp Dp*Cc];%change sign for loop gain
Dout_lqr = -[ Dp*Dc1];

Lout_lqr     = ss(Aout_lqr,Bout_lqr,Cout_lqr,Dout_lqr);
Sout_lqr     = inv(eye(size(Lout_lqr))+Lout_lqr);
Tout_lqr     = eye(size(Lout_lqr))-Sout_lqr;

L_lqr_u      = freqresp(Lin_lqr,w);
L_lqr_mad_dB = 20*log10(squeeze(abs(L_lqr_u)));
wc_lqr           = crosst(L_lqr_mad_dB,w);

RD_lqr_u     = sigma(Lin_lqr,w,2);
SR_lqr_u     = sigma(Lin_lqr,w,3);
rdmin_lqr        = min(RD_lqr_u);
srmin_lqr        = min(SR_lqr_u);

L_lqr_y      = freqresp(Lout_lqr,w);
S_lqr_y      = freqresp(Sout_lqr,w);
T_lqr_y      = freqresp(Tout_lqr,w);
Smax_lqr         = max(abs(S_lqr_y(1,1,:)));
Tmax_lqr         = max(abs(T_lqr_y(1,1,:)));

neg_gm =  min([ (1/(1+rdmin_lqr)) (1-srmin_lqr)]); % in dB
pos_gm =  max([ (1/(1-rdmin_lqr)) (1+srmin_lqr)]); % in dB
neg_gmdB = 20*log10( neg_gm ); % in dB
pos_gmdB = 20*log10( pos_gm );
pm = 180*(max([2*asin(rdmin_lqr/2) 2*asin(srmin_lqr/2)]))/pi;

disp('********Frequency Analysis on lqr************')
LGCF_lqr  =  ['LGCF = ' num2str(wc_lqr) ' rps']
RDMIN_lqr =  ['min|I+L| = ' num2str(rdmin_lqr)]
SRMIN_lqr =  ['min|I+invL| = ' num2str(srmin_lqr)]
SMAX_lqr  =  ['max|S| = ' num2str(Smax_lqr)]
TMAX_lqr  =  ['max|T| = ' num2str(Tmax_lqr)]
disp(['Singular value gain margins = [' ...
         num2str(neg_gmdB) ' dB,' num2str(pos_gmdB) ' dB ]' ])
     disp(['Singular value phase margins = [ +/-' ...
         num2str(pm)  ' deg ]' ])

% connect the controller with plant, form a closed loop system
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
sys_cl = ss(Acl,Bcl,Ccl,Dcl);

% Closed loop eigenvalues and eigenvectors

[Acl_V, Acl_D] = eig(Acl);

%step response analysis
y = step(sys_cl,t);
s = stepinfo(sys_cl);
Az_lqr      = y(:,1);
q_lqr       = y(:,2).*(180/pi);
del_lqr     = (180/pi).*y(:,3);
deldot_lqr  = (180/pi).*y(:,4);

%static projective controller
F = Ar-Br*Kc;
[xx,yy]=eig(F);
lamda = diag(yy);

xr = [ xx(:,3:5)];
xr1 = xr(1:3,1:3);
xr2 = xr(4:5,1:3);
xp = [ real(xx(:,1)) imag(xx(:,1))];
xp1 = xp(1:3,1:2);
xp2 = xp(4:5,1:2);

C= [ eye(3) 0.*ones(3,2)];
ky = real(Kc*xr*inv(C*xr));

%Form the static output feedback controller

Ac_of =  0.;
Bc1_of = [1. 0.];
Bc2_of =  -1;
Cc_of  = -ky(1);
Dc1_of = -ky(2:3);
Dc2_of = 0.;

% plant dimension changed
Ap = [Za/V, Za, 0, Zs;Ma/Za, 0, Ms-(Ma*Zs)/Za, 0;0, 0, 0, 1;0, 0, -wn^2, -2*xi*wn];
Bp = [0; 0; 0; wn^2];
Cp = eye(2,4);
Dp = 0.*Cp*Bp;


Z_of = inv(eye(size(Dc1_of*Dp))-Dc1_of*Dp);
Acl_of = [     (Ap+Bp*Z_of*Dc1_of*Cp)              (Bp*Z_of*Cc_of);
    (Bc1_of*(Cp+Dp*Z_of*Dc1_of*Cp))  (Ac+Bc1_of*Dp*Z_of*Cc_of)];
Bcl_of = [       Bp*Z_of*Dc2_of;
    (Bc2_of+Bc1_of*Dp*Z_of*Dc2_of)];
Ccl_of = [(Cp+Dp*Z_of*Dc1_of*Cp) (Dp*Z_of*Cc_of)];
Dcl_of =(Dp*Z_of*Dc2_of);
sys_clof = ss(Acl_of,Bcl_of,Ccl_of,Dcl_of);

[D,DD] = eig(Acl_of);

%step response 
[y_spc, x_spc] = step(Acl_of,Bcl_of,Ccl_of,Dcl_of,1 ,t);
Az_spc      = y_spc(:,1);
q_spc       = y_spc(:,2).*(180/pi);
del_spc     = x_spc(:,3).*(180/pi);
deldot_spc  = x_spc(:,4).*(180/pi);
int_err_spc = x_spc(:,5);

%frequency analysis of SPC
Ain_SPC = [ Ap 0.*Bp*Cc;  Bc1_of*Cp Ac_of];
Bin_SPC = [ Bp; Bc1_of*Dp];
Cin_SPC = -[ Dc1_of*Cp Cc_of];%change sign for loop gain
Din_SPC = -[ Dc1_of*Dp];

Lin_SPC = ss(Ain_SPC,Bin_SPC,Cin_SPC,Din_SPC);

%SS model of loop gain at the plant output
Aout_SPC = [ Ap Bp*Cc_of;  0.*Bc1_of*Cp Ac_of];
Bout_SPC = [ Bp*Dc1_of; Bc1_of];
Cout_SPC = -[ Cp Dp*Cc_of];%change sign for loop gain
Dout_SPC = -[ Dp*Dc1_of];

Lout_SPC     = ss(Aout_SPC,Bout_SPC,Cout_SPC,Dout_SPC);
Sout_SPC     = inv(eye(size(Lout_SPC))+Lout_SPC);
Tout_SPC     = eye(size(Lout_SPC))-Sout_SPC;

L_spc_u      = freqresp(Lin_SPC,w);
L_spc_mad_dB = 20*log10(squeeze(abs(L_spc_u)));
wc_spc           = crosst(L_spc_mad_dB,w);

RD_spc_u     = sigma(Lin_SPC,w,2);
SR_spc_u     = sigma(Lin_SPC,w,3);
rdmin_spc        = min(RD_spc_u);
srmin_spc        = min(SR_spc_u);

L_spc_y      = freqresp(Lout_SPC,w);
S_spc_y      = freqresp(Sout_SPC,w);
T_spc_y      = freqresp(Tout_SPC,w);
Smax_spc         = max(abs(S_spc_y(1,1,:)));
Tmax_spc         = max(abs(T_spc_y(1,1,:)));

neg_gm_spc =  min([ (1/(1+rdmin_spc)) (1-srmin_spc)]); % in dB
pos_gm_spc =  max([ (1/(1-rdmin_spc)) (1+srmin_spc)]); % in dB
neg_gmdB_spc = 20*log10( neg_gm_spc ); % in dB
pos_gmdB_spc = 20*log10( pos_gm_spc );
pm = 180*(max([2*asin(rdmin_spc/2) 2*asin(srmin_spc/2)]))/pi;

disp('********Frequency Analysis on SPC************')
LGCF_spc  =  ['LGCF = ' num2str(wc_spc) ' rps']
RDMIN_spc =  ['min|I+L| = ' num2str(rdmin_spc)]
SRMIN_spc =  ['min|I+invL| = ' num2str(srmin_spc)]
SMAX_spc  =  ['max|S| = ' num2str(Smax_spc)]
TMAX_spc  =  ['max|T| = ' num2str(Tmax_spc)]

disp(['Singular value gain margins = [' ...
         num2str(neg_gmdB_spc) ' dB,' num2str(pos_gmdB_spc) ' dB ]' ])
     disp(['Singular value phase margins = [ +/-' ...
         num2str(pm)  ' deg ]' ])

%run homework3
RSLQR

%plot the step respose and compare with RSLQR

figure; 
plot(t,Az_lqr,'b','LineWidth',2);  
title('Step Response-Acceleration: 5-state RSLQR V.S. 3-state RSLQR');
ylabel('Az');
xlabel('time(sec)');
grid on
hold on
plot(t,Az1,'g','LineWidth',2);
hold on
% plot(t,Az_spc,'g','LineWidth',2);
legend('5-state RSLQR:rising time = 0.2068','3-state RSLQR:rising time = 0.2068')

%plot the step response of spc and rslqr
figure;
plot(t,Az_lqr,'r','LineWidth',2);
hold on 
grid on
plot(t,Az1,'g','LineWidth',2);
hold on
plot(t,Az_spc,'b','LineWidth',2);
title('Step Response-Acceleration');
legend('5-RSLQR','3-state RSLQR','SPC');
ylabel('Az');
xlabel('time(sec)');

figure;
plot(t,q_lqr,'r','LineWidth',2);
hold on
plot(t,q1,'g','LineWidth',2);
grid on
hold on
plot(t,q_spc,'b','LineWidth',2);
legend('5-RSLQR','3-state RSLQR','SPC');
title('Step Response-pitch rate');
ylabel('q');
xlabel('time(sec)');

figure;
plot(t,del_lqr,'r','LineWidth',2);
hold on
plot(t,del1,'g','LineWidth',2);
hold on
grid on
plot(t,del_spc,'b','LineWidth',2);
legend('5-RSLQR','3-state RSLQR','SPC');
title('Step Response-Elevon');
ylabel('Elevon');
xlabel('time(sec)');

figure;
plot(t,deldot_lqr,'r','LineWidth',2);
hold on
plot(t,deldot1,'g','LineWidth',2);
hold on
grid on
plot(t,deldot_spc,'b','LineWidth',2);
legend('5-RSLQR','3-state RSLQR','SPC');
title('Step Response-Elevon Rate');
ylabel('Elevon Rate');
xlabel('time(sec)');

figure
    semilogx(w,L_lqr_mad_dB,'k',w,L_spc_mad_dB,'b',w,magdb1,'r','LineWidth',2);grid
    title('Bode in dB');
    ylabel('|L| (dB)');
    xlabel('Frequency (rps)');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR');
   
figure
    semilogx(w,rtd*phase(squeeze(L_lqr_u)),'k',w,rtd*phase(squeeze(L_spc_u)),'b',w,rtd*phase(squeeze(Li1)),'r','LineWidth',2);grid
    title('Bode in deg');
    ylabel('Phase (deg)');
    xlabel('Frequency (rps)');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR');
 
figure
    semilogx(w,20*log10(RD_lqr_u),'k',w,20*log10(RD_spc_u),'b',w,20*log10(rd1),'r','LineWidth',2);grid
    title('Return Difference');
    ylabel('|I+L| (dB)');
    xlabel('Frequency (rps)');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR' );
 
figure
    semilogx(w,20*log10(SR_lqr_u),'k',w,20*log10(SR_spc_u),'b',w,20*log10(sr1),'r','LineWidth',2);grid
    title('Stability Robustness');
    ylabel('|I+invL| (dB)');
    xlabel('Frequency (rps)');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR' );
 
figure
    th=0.:.1:2*pi;
    plot(real(squeeze(L_lqr_u)),imag(squeeze(L_lqr_u)),'k',...
         real(squeeze(L_spc_u)),imag(squeeze(L_spc_u)),'b',...
         real(squeeze(Li1)),imag(squeeze(Li1)),'g',... 
           sin(th),cos(th),'r','LineWidth',2);grid
    axis([-2 2 -2 2]);
    title('Nyquist');
    xlabel('Real');ylabel('Imag');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR','Unit Circle' );

figure
    semilogx(w,20*log10(squeeze(abs(S_lqr_y(1,1,:)))),'k',w,20*log10(squeeze(abs(S_spc_y(1,1,:)))),'b',w,20*log10(squeeze(abs(Slo_y(1,1,:)))),'r','LineWidth',2);grid
    title('Sensitivity Az ');
    ylabel('|S| (dB)');
    xlabel('Frequency (rps)');
    legend('5-State RSLQR' ,'SPC','3-State RSLQR');
  
figure
    semilogx(w,20*log10(squeeze(abs(T_lqr_y(1,1,:)))),'k',w,20*log10(squeeze(abs(T_spc_y(1,1,:)))),'b',w,20*log10(squeeze(abs(Tlo_y(1,1,:)))),'r','LineWidth',2);grid
    title('Complementary Sensitivity Az');
    ylabel('|T| (dB)');
    xlabel('Frequency (rps)');

legend('5-State RSLQR' ,'SPC','3-State RSLQR' );

%********Matrices List*********************
disp('***************Awiggle is:')
Ar
disp('***************Bwiggle is:')
Br
disp('***************Q_lqr is:')
Q
disp('***************R_lqr is:')
R
disp('***************Kx_lqr is:')
Kc
disp('***************Controller Matrices are as follows:')
Ac
Bc1
Bc2
Cc
Dc1
Dc2
disp('***********Closed loop eigenvalues in diagonal matrix and eigenvectors(each column)************')
Acl_V
Acl_D

disp('***********Static Projetive Controller')
F
disp('***************Output Feedback Controller Matrices are as follows:')
Ac_of
Bc1_of
Bc2_of
Cc_of
Dc1_of
Dc2_of
disp('***********Eigenvalue and Eigenvectors of Output feedback ')
D
DD





