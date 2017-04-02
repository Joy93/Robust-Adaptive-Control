%% Homework4
% Yejing Zhou & Mengyan Li

w = logspace(-3,4,1000);
t = linspace(0,2.5,1000);
%% Plant Model
Ap = [-1.3046e+00   1.0000e+00  -2.1420e-01            0;
   4.7711e+01            0  -1.0483e+02            0;
            0            0            0   1.0000e+00;
            0            0  -4.7769e+03  -9.7729e+01]


Bp =[0;0;0;4.7769e+03]


Cp =[ -1.1569e+03            0  -1.8995e+02            0;
   1.0000e+00            0            0            0;
            0   1.0000e+00            0            0;
            0            0   1.0000e+00            0;
            0            0            0   1.0000e+00]


Dp =[0;
     0;
     0;
     0;
     0]

%% Target Loop Gain Cross-Over Frequency
% This parameter set the target bandwidth

wlgcf = 1.2; %Hz
disp(' ')
disp(['Target lgcf = ' num2str(wlgcf)])
disp(' ')

%% Sensitivity Ws

Ws.Tau = 1/(2*pi*wlgcf);
Ws.K = 0.5/Ws.Tau;
%Ws.K = 0.8/Ws.Tau;
Ws.num = Ws.K * [Ws.Tau, 1];
Ws.den = [1, 0];
[Ws.A,Ws.B,Ws.C,Ws.D] = tf2ss(Ws.num,Ws.den);

Ws.sys = ss(Ws.A,Ws.B,Ws.C,Ws.D);
figure(1)
bodemag(Ws.sys)
grid on
hold on

%% Complementary Sensitivity Wt

Wt.TauN = 1/(3.5*2*pi*wlgcf);
Wt.TauD = 0.005;
Wt.K = 0.707;
Wt.num = Wt.K * [Wt.TauN, 1];
Wt.den = [Wt.TauD, 1];
[Wt.A,Wt.B,Wt.C,Wt.D] = tf2ss(Wt.num,Wt.den);

Wt.sys = ss(Wt.A,Wt.B,Wt.C,Wt.D);


%% Control Activity

Wc = 0.1;
Cc = [0, 0, 0, 1];
Wc_sys = ss(0,0,0,Wc);

%% Controller
%%Hinf Controller
Ac_hinf = [0     0; 0  -200]
Bc1_hinf = [ 1     0     0     0     0; 1     0     0     0     0]
Bc2_hinf = [-1;0]
Cc_hinf =[-2.0470e-02  -1.4680e-01]
Dc1_hinf = [0   5.1027e+00   4.6076e-01  -7.0538e-01  -6.3309e-03]
Dc2_hinf = 0

%RSLQR Controller
Ac_r = 0
Bc1_r = [1     0     0     0     0]
Bc2_r = -1
Cc_r = -9.3194e-03
Dc1_r = [0   2.5923e+00   2.3625e-01   0            0]
Dc2_r = 0;

%% Closed loop system
%Closed loop System for Hinf
 Z_hinf = inv(eye(size(Dc1_hinf*Dp))-Dc1_hinf*Dp);
 Acl_hinf = [(Ap+Bp*Z_hinf*Dc1_hinf*Cp)              (Bp*Z_hinf*Cc_hinf);
         (Bc1_hinf*(Cp+Dp*Z_hinf*Dc1_hinf*Cp))  (Ac_hinf+Bc1_hinf*Dp*Z_hinf*Cc_hinf)]
 Bcl_hinf = [       Bp*Z_hinf*Dc2_hinf;
         (Bc2_hinf+Bc1_hinf*Dp*Z_hinf*Dc2_hinf)]
 Ccl_hinf = [(Cp+Dp*Z_hinf*Dc1_hinf*Cp) (Dp*Z_hinf*Cc_hinf)]
 Dcl_hinf=(Dp*Z_hinf*Dc2_hinf)
 sys_clh = ss(Acl_hinf,Bcl_hinf,Ccl_hinf,Dcl_hinf);

 %Closed loop system for RSLQR
 Z_r = inv(eye(size(Dc1_r*Dp))-Dc1_r*Dp);
 Acl_r = [     (Ap+Bp*Z_r*Dc1_r*Cp)              (Bp*Z_r*Cc_r);
         (Bc1_r*(Cp+Dp*Z_r*Dc1_r*Cp))  (Ac_r+Bc1_r*Dp*Z_r*Cc_r)]
 Bcl_r = [       Bp*Z_r*Dc2_r;
         (Bc2_r+Bc1_r*Dp*Z_r*Dc2_r)]
 Ccl_r = [(Cp+Dp*Z_r*Dc1_r*Cp) (Dp*Z_r*Cc_r)]
 Dcl_r =(Dp*Z_r*Dc2_r)
 sys_clr = ss(Acl_r,Bcl_r,Ccl_r,Dcl_r);

%% Step Response
 %Step response for Hinf
 y_hinf = step(sys_clh,t);
 az_h = y_hinf(:,1); %  acceleration (fps2)
 aze_h = abs(ones(size(az_h))-az_h);  % error for az
 taur_h = 0.; taus_h= 0.; % rise time and settling time
 fv_h = aze_h(numel(aze_h)); % final value of the error
 e_nh = aze_h - fv_h*ones(size(aze_h)) - 0.36*ones(size(aze_h));
 e_n1h = abs(e_nh) + e_nh;
 taur_h = crosst(e_n1h,t); % rise time 
 e_nh = aze_h - fv_h*ones(size(aze_h)) - 0.05*ones(size(aze_h));
 e_n1h = abs(e_nh) + e_nh;
 taus_h = crosst(e_n1h,t); % settling time
%Step response for RSLQR
  y_r = step(sys_clr,t);
  az_r = y_r(:,1); %  acceleration (fps2)
  aze_r = abs(ones(size(az_r))-az_r);  % error for az
  taur_r = 0.; taus_r= 0.; % rise time and settling time
  fv_r = aze_r(numel(aze_r)); % final value of the error
  e_nr = aze_r - fv_r*ones(size(aze_r)) - 0.36*ones(size(aze_r));
  e_n1r = abs(e_nr) + e_nr;
  taur_r = crosst(e_n1r,t); % rise time 
  e_nr = aze_r - fv_r*ones(size(aze_r)) - 0.05*ones(size(aze_r));
  e_n1r = abs(e_nr) + e_nr;
  taus_r = crosst(e_n1r,t); % settling time
     

%% Frequency Domain Analysis
    %RSLQR
  A_Lur = [ Ap 0.*Bp*Cc_r;  Bc1_r*Cp Ac_r];
  B_Lur = [ Bp; Bc1_r*Dp];
  C_Lur = -[ Dc1_r*Cp Cc_r];%change sign for loop gain
  D_Lur = -[ Dc1_r*Dp];
  sys_Lur = ss(A_Lur,B_Lur,C_Lur,D_Lur);
  magdb = 20*log10(abs(squeeze(freqresp(sys_Lur,w))));
  wc_r = crosst(magdb,w); % LGCF, assumes Lu is a scalar
  sr_r = sigma(sys_Lur,w,3);
  sru_minr = min(abs(sr_r));
  rd_r = sigma(sys_Lur,w,2);
  rdu_minr = min(abs(rd_r));
  Lu_r = freqresp(sys_Lur,w);
    %Hinf
  A_Luh = [ Ap 0.*Bp*Cc_hinf;  Bc1_hinf*Cp Ac_hinf];
  B_Luh = [ Bp; Bc1_hinf*Dp];
  C_Luh = -[ Dc1_hinf*Cp Cc_hinf];%change sign for loop gain
  D_Luh = -[ Dc1_hinf*Dp];
  sys_Luh = ss(A_Luh,B_Luh,C_Luh,D_Luh);
  magdb = 20*log10(abs(squeeze(freqresp(sys_Luh,w))));
  wc_h = crosst(magdb,w); % LGCF, assumes Lu is a scalar
  sr_h = sigma(sys_Luh,w,3);
  sru_minh = min(abs(sr_h));
  rd_h = sigma(sys_Luh,w,2);
  rdu_minh = min(abs(rd_h));
  Lu_h = freqresp(sys_Luh,w);
     
%% Plots
% Plot Unit step Az command
figure;
plot(t,az_h,'b','LineWidth',2);
hold on
plot(t,az_r,'r','LineWidth',2);
grid
title('Step Reponse of Acceleration: Hinf V.S. RSLQR')
legend(['Hinf-Response: 63% Tr = ' num2str(taur_h) ' 95% Ts = ' num2str(taus_h)],['RSLQR-Response: 63% Tr = ' num2str(taur_r) ' 95% Ts = ' num2str(taus_r)])
xlabel('Time (sec)');
ylabel('Az (fps2)');

% Plot Nyquist
  theta = -2*pi:.001:2*pi;
  x1 = cos(theta);
  y1 = sin(theta);
  xu = x1-1;
  yu = y1;
  xir = rdu_minr*x1 - 1;
  yir = rdu_minr*y1 ;
  xih = rdu_minh*x1 - 1;
  yih = rdu_minh*y1 ;
  figure 
  plot(xih,yih,'-.r',xir,yir,'r',xu,yu,'k',real(squeeze(Lu_r)),imag(squeeze(Lu_r)),'b',real(squeeze(Lu_h)),imag(squeeze(Lu_h)),'g');
  grid
  axis([-3 2 -2 2]);
  legend('Loop-Hinf','Loop-RSLQR','Unit Circle','Lu-RSLQR','Lu-Hinf');
  xlabel('Real')
  ylabel('Imaginary')
  title('Nyquist Plot at Plant Inutput: RSLQR V.S. Hinf')

%plot bode
  figure;
  margin(sys_Lur);grid
  hold on
  margin(sys_Luh)
  legend('RSLQR:','Hinf')
  xlabel('Frequency (rps)')
  ylabel('Mag')
  title('Bode at Input: RSLQR V.S. Hinf')
%plot return difference
figure
semilogx(w,20*log10(abs(rd_r)),'r','LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(rd_h)),'b','LineWidth',2);
legend([' RSLQR-min(I+Lu_r) = ' num2str(rdu_minr)],[' Hinf-min(I+Lu_h) = ' num2str(rdu_minh)]);
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input: RSLQR V.S. Hinf')
%plot stability matrix
figure
semilogx(w,20*log10(abs(sr_r)),'r','LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(sr_h)),'b','LineWidth',2);
legend([' RSLQR-min(I+Lu_r) = ' num2str(sru_minr)],[' Hinf-min(I+Lu_h) = ' num2str(sru_minh)]);
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input: RSLQR V.S. Hinf')
  
%%Compute singluar value margins
%H-inf
hinf_neg_gm =  min([ (1/(1+rdu_minh)) (1-sru_minh)]); % in dB
hinf_pos_gm =  max([ (1/(1-rdu_minh)) (1+sru_minh)]); % in dB
hinf_neg_gmdB = 20*log10( hinf_neg_gm ); % in dB
hinf_pos_gmdB = 20*log10( hinf_pos_gm ); % in dB
hinf_pm = 180*(max([2*asin(rdu_minh/2) 2*asin(sru_minh/2)]))/pi;% in deg
disp('H inf Singular value margins')
disp(['H inf Min Singular value I+Lu =    ' num2str(rdu_minh)])
disp(['H inf Min Singular value I+invLu = ' num2str(sru_minh)])
disp(['H inf Singular value gain margins = [' ...
num2str(hinf_neg_gmdB) ' dB,' num2str(hinf_pos_gmdB) ' dB ]' ])
disp(['H inf Singular value phase margins = [ +/-' ...
num2str(hinf_pm)  ' deg ]' ])

%RSLQR
rslqr_neg_gm =  min([ (1/(1+rdu_minr)) (1-sru_minr)]); % in dB
rslqr_pos_gm =  max([ (1/(1-rdu_minr)) (1+sru_minr)]); % in dB
rslqr_neg_gmdB = 20*log10( rslqr_neg_gm ); % in dB
rslqr_pos_gmdB = 20*log10( rslqr_pos_gm ); % in dB
rslqr_pm = 180*(max([2*asin(rdu_minr/2) 2*asin(sru_minr/2)]))/pi;% in deg
disp('RSLQR singular value margins')
disp(['RSLQR Min Singular value I+Lu =    ' num2str(rdu_minr)])
disp(['RSLQR Min Singular value I+invLu = ' num2str(sru_minr)])
disp(['RSLQR Singular value gain margins = [' ...
num2str(rslqr_neg_gmdB) ' dB,' num2str(rslqr_pos_gmdB) ' dB ]' ])
disp(['RSLQRSingular value phase margins = [ +/-' ...
num2str(rslqr_pm)  ' deg ]' ])
     
% Analysis at plant output
%H-inf     
hinf_T  = freqresp(sys_clh,w); % Complementary Sensitivity
hinf_S = 1 - hinf_T; % Sensitivity
hinf_T_Az = 20*log10(abs(squeeze(hinf_T(1,1,:))));
hinf_S_Az = 20*log10(abs(squeeze(hinf_S(1,1,:))));
hinf_Tmax = max(hinf_T_Az); % Inf Norm of T in dB
hinf_Smax = max(hinf_S_Az); % Inf Norm of S in dB
hinf_Ws_magdb = 20*log10(abs(squeeze(1/freqresp(Ws.sys,w))));
hinf_Wt_magdb = 20*log10(abs(squeeze(1/freqresp(Wt.sys,w)))); 

%RSLQR
rslqr_T  = freqresp(sys_clr,w); % Complementary Sensitivity
rslqr_S = 1 - rslqr_T; % Sensitivity
rslqr_T_Az = 20*log10(abs(squeeze(rslqr_T(1,1,:))));
rslqr_S_Az = 20*log10(abs(squeeze(rslqr_S(1,1,:))));
rslqr_Tmax = max(rslqr_T_Az); % Inf Norm of T in dB
rslqr_Smax = max(rslqr_S_Az); % Inf Norm of S in dB
% rslqr_Ws_magdb = 20*log10(abs(squeeze(1/freqresp(Ws.sys,w)))) ;
% rslqr_Wt_magdb = 20*log10(abs(squeeze(1/freqresp(Wt.sys,w)))) ;

figure('Name','Comp Sens T');
semilogx(w,hinf_T_Az,'b', w,rslqr_T_Az,'r', w,hinf_Wt_magdb,'k','LineWidth',2);grid
legend(['H-inf: ||T||inf = ' num2str(hinf_Tmax) ' (dB)'],['RSLQR: ||T||inf = ' num2str(rslqr_Tmax) ' (dB)'], 'invWt','Location','Best');
title('Comp Sens T');
xlabel('Freq (rps)');ylabel('Mag (dB)');


figure('Name','Sens S');
semilogx(w,hinf_S_Az,'b',w,rslqr_S_Az,'r',w,hinf_Ws_magdb,'k','LineWidth',2);grid
legend(['H-inf: ||S||inf = ' num2str(hinf_Smax) ' (dB)'],['RSLQR: ||S||inf = ' num2str(rslqr_Smax) ' (dB)'],'invWs','Location','Best');
title('Sens S');
xlabel('Freq (rps)');ylabel('Mag (dB)');
