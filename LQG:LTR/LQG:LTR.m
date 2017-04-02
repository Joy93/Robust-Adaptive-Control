%Homework6
%Yejing Zhou && Mengyan Li
%**************************************************************
clc
clear
% set plot scales
rtd  = 180/pi;
t = 0:0.002:5;
w = logspace(-1,3,500);
dd = 0:0.001:2*pi;
xx1 = cos(dd) - 1;
yy1 = sin(dd);
% actuator model
wn = 2*pi*11;
xi = 0.707;
%plant model
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
% plant model for design without actuator
Ad = [Za/V 1.; Ma 0.];
Bd = [Zs/V;Ms];

% Define wiggle matrices
Cw = [1. 0. 0.; 0. 0. 1.];
Aw = [0, Za, 0; 0, Za/V, 1; 0, Ma, 0];
Bw = [Zs;Zs/V;Ms];

Gw = eye(3);

% Solve Klqr
q=0.8685e-4*[  1    0       0   ;
       0    0       0    ;
       0    0       0       ];
[Kc,~,~]=lqr(Aw,Bw,q,1.);
disp('***************************Gain Matrix for LQR*********************')
Kc
% LQR Analysis
Ac = 0.;
Bc1 = [1. 0. 0. 0. 0.];
Bc2 = -1;
Cc = -Kc(1);
Dc1 = [0. -Kc(2:3) 0. 0.];
Dc2 = 0.;

% Closed-Loop LQR
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);
(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)]; 
Bcl = [ Bp*Z*Dc2;
(Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)]; 
Dcl =(Dp*Z*Dc2);
sys_cl_lqr = ss(Acl,Bcl,Ccl,Dcl);

% RSLQR Analysis
disp('*****************Sigular Value for LQR****************************')
A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
B_Lu = [ Bp; Bc1*Dp];
C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
D_Lu = -[ Dc1*Dp];
sys_Lu_lqr = ss(A_Lu,B_Lu,C_Lu,D_Lu);
Lu_lqr = freqresp(sys_Lu_lqr,w);
magdb_lqr = 20*log10(abs(squeeze(Lu_lqr)));
a = angle(squeeze(Lu_lqr));
phs_lqr = a - (sign(a)==1)*2*pi;
wc_lqr = crosst(magdb_lqr,w); % LGCF, assumes Lu is a scalar
wg_lqr = crosst(phs_lqr+pi,w); %PCF
sr_lqr = sigma(sys_Lu_lqr,w,3);
sru_min_lqr = min(abs(sr_lqr))
rd_lqr = sigma(sys_Lu_lqr,w,2);
rdu_min_lqr = min(abs(rd_lqr))

% Analysis at plant output
T  = freqresp(sys_cl_lqr,w); % Complementary Sensitivity
S = 1 - T; % Sensitivity
T_Az_lqr = 20*log10(abs(squeeze(T(1,1,:))));
S_Az_lqr = 20*log10(abs(squeeze(S(1,1,:))));
Tmax_lqr = max(T_Az_lqr); % Inf Norm of T in dB
Smax_lqr = max(S_Az_lqr); % Inf Norm of S in dB

%SS model of loop gain at the plant output
Av = [Ap+Bp*Z*Dc1*Cp                          Bp*Z*Cc;
      Bc1*(eye(size(Dp*Z*Dc1))+Dp*Z*Dc1)*Cp   Ac+Bc1*Dp*Z*Cc];
Bv = [Bp*Z*Dc1;   Bc1+Bc1*Dp*Z*Dc1];
Cv = [Z*Dc1*Cp   Z*Cc];
Dv = Z*Dc1;

for a=1:numel(w),
    s = sqrt(-1)*w(a);
    KK = Cv*inv(s*eye(size(Av))-Av)*Bv+Dv;
    KKr = s*KK;
    Gnois(a)   = max(svd(KK));
    Gnoisr(a) = max(svd(KKr));
end

glqr = Gnois;
grlqr = Gnoisr;

%Compute singluar value margins
disp('******************Gain Margin and Phase Margin of LQR**************')
neg_gm =  min([ (1/(1+rdu_min_lqr)) (1-sru_min_lqr)]) % in dB
pos_gm =  max([ (1/(1-rdu_min_lqr)) (1+sru_min_lqr)]) % in dB
neg_gmdB = 20*log10( neg_gm ) % in dB
pos_gmdB = 20*log10( pos_gm ) % in dB
pm = 180*(max([2*asin(rdu_min_lqr/2) 2*asin(sru_min_lqr/2)]))/pi % in deg
ngm_lqr = -1/neg_gm
pgm_lqr = -1/pos_gm



% Time Domain LQR
y_lqr = step(sys_cl_lqr,t);
az = y_lqr(:,1); %  Az (fps2)
aze = abs(ones(size(az))-az);  % error for Az
taur_lqr = 0.; % rise time
taus_lqr= 0.; %settling time
fv = aze(numel(aze)); % final value of the error
e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taur_lqr = crosst(e_n1,t); % rise time 
e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taus_lqr = crosst(e_n1,t); % settling time

%Q and R matrix
Q0 = diag([ 0.001 0.0014 0.005]);
R0 = 1000*diag([ 0.025^2 0.001^2]);
p = [10^10 10^4 10^2];
for i = 1:1:size(p,2) 
    Q = Q0+Bw*Bw'/p(i);
    R = R0;
    [Lp,Pf] = lqe(Aw,Gw,Cw,Q,R);
    Ac =  [0.        zeros(1,3);
           Lp(:,1)   Aw-Bw*Kc-Lp*Cw]; 
    Bc1 = [1.   0.      0.        0.   0.;
           zeros(3,2)   Lp(:,2)   zeros(3,2);];
    Bc2 = [-1;   -1;   0.;   0.;];
    Cc = [0 -Kc];
    Dc1 = zeros(1,5);
    Dc2 = 0.;
    
% closed-loop system
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
sys_cl = ss(Acl, Bcl, Ccl, Dcl);

%SS model of loop gain at the plant input
Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc];%change sign for loop gain
Din = -[ Dc1*Dp];
L_input = ss(Ain,Bin,Cin,Din);

% Time Domain
y = step(sys_cl,t);
az = y(:,1); %  Az (fps2)
aze = abs(ones(size(az))-az);  % error for Az
fv = aze(numel(aze)); % final value of the error
e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taur(i) = crosst(e_n1,t); % rise time 
e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
e_n1 = abs(e_n) + e_n;
taus(i) = crosst(e_n1,t); % settling time    
yp(:,:,i)=y;

disp('******************Singular Value for LQG**************************')
A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
B_Lu = [ Bp; Bc1*Dp];
C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
D_Lu = -[ Dc1*Dp];
sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
Lu = freqresp(sys_Lu,w);
magdb = 20*log10(abs(squeeze(Lu)));
a = angle(squeeze(Lu));
phs = a - (sign(a)==1)*2*pi;
wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
wg = crosst(phs+pi,w); %PCF
sr = sigma(sys_Lu,w,3);
sru_min = min(abs(sr))
rd = sigma(sys_Lu,w,2);
rdu_min = min(abs(rd))
 
magdbp(:,i) = magdb;
phsp(:,i) = phs;
sLu(:,i) = squeeze(Lu);
sys_Lup(:,i) = sys_Lu;
srp(:,i) = sr;
sru(i) = sru_min;
rdp(:,i) = rd;
rdu(i) = rdu_min;
    
% Analysis at plant output
T  = freqresp(sys_cl,w); % Complementary Sensitivity
S = 1 - T; % Sensitivity
T_Az = 20*log10(abs(squeeze(T(1,1,:))));
S_Az = 20*log10(abs(squeeze(S(1,1,:))));
Tmax = max(T_Az); % Inf Norm of T in dB
Smax = max(S_Az); % Inf Norm of S in dB    
T_Azp(:,i) = T_Az;
S_Azp(:,i) = S_Az;
Tmaxp(i) = Tmax;
Smaxp(i) = Smax;

%Noise
Av = [Ap+Bp*Z*Dc1*Cp                          Bp*Z*Cc;
      Bc1*(eye(size(Dp*Z*Dc1))+Dp*Z*Dc1)*Cp   Ac+Bc1*Dp*Z*Cc];
Bv = [Bp*Z*Dc1;   Bc1+Bc1*Dp*Z*Dc1];
Cv = [Z*Dc1*Cp   Z*Cc];
Dv = Z*Dc1;

for a=1:numel(w)
    s = sqrt(-1)*w(a);
    KK = Cv*inv(s*eye(size(Av))-Av)*Bv+Dv;
    KKr = s*KK;
    Gnois(a)   = max(svd(KK));
    Gnoisr(a) = max(svd(KKr));
end

LGn(:,i) = Gnois;
LGnr(:,i) = Gnoisr;

%Compute singluar value margins
disp('******************Gain Margin and Phase Margin of LQG**************')
neg_gm =  min([ (1/(1+rdu_min)) (1-sru_min)]) % in dB
pos_gm =  max([ (1/(1-rdu_min)) (1+sru_min)]) % in dB
neg_gmdB = 20*log10( neg_gm ) % in dB
pos_gmdB = 20*log10( pos_gm ) % in dB
pm = 180*(max([2*asin(rdu_min/2) 2*asin(sru_min/2)]))/pi% in deg
ngm = -1/neg_gm
pgm = -1/pos_gm

end
% Wiggle Matrices
disp('**********************Observer Design Matrices********************')
Aw
Bw
Cw
Gw
Q
R
Lp
disp('**********************Controller Design Matrices********************')
Ac
Bc1
Bc2
Cc
Dc1
Dc2




% Time plot
     
figure
plot(t,y_lqr(:,1),t,yp(:,1,1),t,yp(:,1,2),t,yp(:,1,3),'LineWidth',2);grid
legend(['RSLQR 63% Tr = ' num2str(taur_lqr) ' 95% Ts = ' num2str(taus_lqr)],['p = ' num2str(p(1)) ' 63% Tr = ' num2str(taur(1)) ' 95% Ts = ' num2str(taus(1))],['p = ' num2str(p(2)) ' 63% Tr = ' num2str(taur(2)) ' 95% Ts = ' num2str(taus(2))],['p = ' num2str(p(3)) ' 63% Tr = ' num2str(taur(3)) ' 95% Ts = ' num2str(taus(3))]);
xlabel('Time (sec)');
ylabel('Az (fps2)');
title('Step Response-Az');

figure
plot(t,y_lqr(:,2),t,yp(:,2,1),t,yp(:,2,2),t,yp(:,2,3),'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Time (sec)');
ylabel('AOA (rad)');
title('Step Response-AOA');

figure('Name','Pitch Rate Time History')
plot(t,y_lqr(:,3)*rtd,t,yp(:,3,1)*rtd,t,yp(:,3,2)*rtd,t,yp(:,3,3)*rtd,'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Time (sec)');
ylabel('Pitch Rate (dps)');
title('Step Response-Pitch Rate');

figure
plot(t,y_lqr(:,4)*rtd,t,yp(:,4,1)*rtd,t,yp(:,4,2)*rtd,t,yp(:,4,3)*rtd,'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Time (sec)');
ylabel('Elevon (deg)');
title('Step Response-Elevon');

figure
plot(t,y_lqr(:,5)*rtd,t,yp(:,5,1)*rtd,t,yp(:,5,2)*rtd,t,yp(:,5,3)*rtd,'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Time (sec)');
ylabel('Elevon Rate(dps)');
title('Step Response-Elevon Rate');


% Freq Plot

figure
semilogx(w,magdb_lqr,w,magdbp,'LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title('Bode Magnitude at Plant Input');
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);


figure
semilogx(w,rtd*phs,w,rtd*phsp,'LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Phase deg')
title('Bode Phase at Plant Input');
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);


figure
plot(xx1,yy1,'k:',real(squeeze(Lu_lqr)),imag(squeeze(Lu_lqr)),real(sLu),imag(sLu),...
    [ngm -1.],[0. 0.],'r',[-1 pgm],[0. 0.],'c','LineWidth',2);grid
axis([-3 3 -3 3]);
xlabel('Re(Lu)')
ylabel('Im(Lu)')
legend('Unit Circle at -1,j0','RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))],'Neg SV Margin for p = 10^2','Pos SV Margin for p = 10^2');
title('Nyquist Plot at Plant Input');


figure
nyquist(sys_Lu_lqr,sys_Lup(1),sys_Lup(2),sys_Lup(3))
axis([-2 2 -2 2]);
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot in Plant Input');


figure
semilogx(w,20*log10(abs(rd_lqr)),w,20*log10(abs(rdp)),'LineWidth',2);grid
legend([' RSLQR min(I+Lu) = ' num2str(rdu_min_lqr)],['p = ' num2str(p(1)) ' min(I+Lu) = ' num2str(rdu(1))],['p = ' num2str(p(2)) ' min(I+Lu) = ' num2str(rdu(2))],['p = ' num2str(p(3)) ' min(I+Lu) = ' num2str(rdu(3))],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')


figure
semilogx(w,20*log10(abs(sr_lqr)),w,20*log10(abs(srp)), 'LineWidth',2);grid
legend(['RSLQR min(I+invLu) = ' num2str(sru_min_lqr)],['p = ' num2str(p(1)) ' min(I+invLu) = ' num2str(rdu(1))],['p = ' num2str(p(2)) ' min(I+invLu) = ' num2str(rdu(2))],['p = ' num2str(p(3)) ' min(I+invLu) = ' num2str(rdu(3))],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')


figure
semilogx(w,T_Az_lqr,w,T_Azp,'LineWidth',2);grid
legend(['RSLQR ||T||inf = ' num2str(Tmax_lqr) ' (dB)'],['p = ' num2str(p(1)) ' ||T||inf =  ' num2str(Tmaxp(1)) '(dB)' ],['p = ' num2str(p(2)) ' ||T||inf =  ' num2str(Tmaxp(2)) '(dB)' ],['p = ' num2str(p(3)) ' ||T||inf =  ' num2str(Tmaxp(3)) '(dB)' ],'Location','Best');
title('Compensative Sensitivity T');
xlabel('Freq (rps)');
ylabel('Mag (dB)');


figure
semilogx(w,S_Az,w,S_Azp,'LineWidth',2);grid
legend(['RSLQR ||S||inf = ' num2str(Smax_lqr) ' (dB)'],['p = ' num2str(p(1)) ' ||S||inf =  ' num2str(Smaxp(1)) '(dB)' ],['p = ' num2str(p(2)) ' ||S||inf =  ' num2str(Smaxp(2)) '(dB)' ],['p = ' num2str(p(3)) ' ||S||inf =  ' num2str(Smaxp(3)) '(dB)' ],'Location','Best');
title('Sensitivity S');
xlabel('Freq (rps)');
ylabel('Mag (dB)');

figure
semilogx(w,20*log10(glqr),'b',w,20*log10(LGn),'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Frequency (rps)')
ylabel('Magnitude (dB)')
title('Noise-to-Control Gain Matrix')

figure
semilogx(w,20*log10(grlqr),'b',w,20*log10(LGnr),'LineWidth',2);grid
legend('RSLQR',['p = ' num2str(p(1))],['p = ' num2str(p(2))],['p = ' num2str(p(3))]);
xlabel('Frequency (rps)')
ylabel('Magnitude (dB)')
title('Noise-to-Control Rate Gain Matrix')

