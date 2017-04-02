%Homework 4
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
Ap_hw3 = [-1.3046,1.0, -0.21420, 0; 47.711, 0, -104.83, 0;0, 0, 0, 1; 0, 0, -wn^2, -2*xi*wn];
Bp_hw3 = [0;0;0;wn^2];
Cp_hw3 = [-1156.9, 0, -189.95, 0; eye(4)];
Dp = 0.*Cp_hw3*Bp_hw3;
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
for ii = 1:numel(qq),
    Q(1,1) = qq(ii);
    pen(ii) = qq(ii);
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
    Acl = [(Ap_hw3 + Bp_hw3*z*Dc1*Cp_hw3) (Bp_hw3*z*Cc); (Bc1*(Cp_hw3+Dp*z*Dc1*Cp_hw3)) (Ac+Bc1*Dp*z*Cc)];
    Bcl = [Bp_hw3*z*Dc2; (Bc2+Bc1*Dp*z*Dc2)];
    Ccl = [(Cp_hw3+Dp*z*Dc1*Cp_hw3) (Dp*z*Cc)];
    Dcl = (Dp*z*Dc2);
    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
   
    %compute the eigenvalues for Acl to plot root locus plot
    xacl = eig(Acl);
    xeig = [xeig;xacl];
    
    %loop gain model at input breakpoint
    Ali = [Ap_hw3, 0.*Bp_hw3*Cc; Bc1*Cp_hw3, Ac];
    Bli = [Bp_hw3; Bc1*Dp];
    Cli = -[Dc1*Cp_hw3, Cc];
    Dli = -[Dc1*Dp];
    sys_li = ss(Ali,Bli,Cli,Dli);
    %frequency analysis********************
    Li = freqresp(sys_li,w);
    T = freqresp(sys_cl,w);
    S = 1 - T;
    magdb = 20*log10(abs(squeeze(Li)));
    
    sr = sigma(sys_li,w,3);
    sr_min = min(abs(sr)); %1+invLi
    sr_st(ii) = sr_min;
    rd = sigma(sys_li,w,2);
    rd_min = min(abs(rd)); % 1+Li
    rd_st(ii) = rd_min; %when ii = 21, we find sigma(1+Lu) closest to 0.5.
    
    if abs(rd_min-sr_min) <= 1e-5,
        fi = ii
        abs(rd_min - sr_min)
        break
    end
    step(sys_cl)
    y = step(sys_cl,t);
    Az = y(:,1);
    q = y(:,3).*(180/pi);
    Aze = abs(ones(size(Az)) - Az);
    fe = Aze(numel(Aze));%final error value
    
    Az_st(ii,:) = Az';
    q_st(ii,:) = q;
    Del_st(ii,:) = (180/pi)*y(:,4);
    Deldot_st(ii,:) = (180/pi)*y(:,5);
    
%     if ii == 20,
%         display('*********Matrices in Design************')
%         Ap
%         Bp
%         Cp
%         Dp
%         Acl
%         Bcl
%         Ccl
%         Dcl
%         Ac
%         Bc1
%         Bc2
%         Cc
%         Dc1
%         Dc2
%         display('*********Pick sigma(1+Lu) = 0.5 with following parameters*******')
%         Q
%         Kx_lqr
%         S = stepinfo(sys_cl);
%         display('*******Step Response Information for Acceleration********')
%         display(S(1,:));
%     end
%     
%     
%  
end
% %plot the step response
% 
% 
figure; 
for ii = 1:fi,
    plot(t,Az_st(ii,:));%acceleration
    hold on
end
% title('Step Response-Acceleration');
% ylabel('Az');
% xlabel('time(sec)');
% grid on
% figure; 
% for ii = 1:fi,
%     plot(t,q_st(ii,:));%pitch rate
%     hold on
% end
% title('Step Response-Pitch Rate');
% ylabel('q');
% xlabel('time(sec)');
% grid on 
% figure; 
% for ii = 1:fi,
%     plot(t,Del_st(ii,:));%elevon
%     hold on
% end
% title('Step Response-Elevon');
% ylabel('Elevon');
% xlabel('time(sec)');
% grid on
% figure; 
% for ii = 1:fi,
%     plot(t,Deldot_st(ii,:));%elevon rate
%     hold on
%end
% title('Step Response-Elevon Rate');
% ylabel('Elevon Rate');
% xlabel('time(sec)');
% grid on
% 
% 
% 
% %plot root locus
% figure;
% plot(xeig, 'x');
% hold on
% plot(xol,0, 'o');
% grid on;
% title('RSLQR and Close Loop Root Locus')
% axis([-20 10 -20 20]);
% xlabel('Real');
% ylabel('Imaginary');
% 
% %choose ii = 20
% figure;
% suptitle('Unit Step Response for qq = 0.8685e-4 ')
% subplot(2,2,1);
% plot(t,Az_st(20,:));
% xlabel('time(sec)');
% ylabel('Acceleration');
% grid on
% subplot(2,2,2)
% plot(t,q_st(20,:));
% xlabel('time(sec)');
% ylabel('Pitch Rate');
% grid on
% subplot(2,2,3)
% plot(t,Del_st(20,:));
% xlabel('time(sec)');
% ylabel('Elevon');
% grid on
% subplot(2,2,4)
% plot(t,Deldot_st(20,:));
% xlabel('time(sec)');
% ylabel('Elevon Rate');
% grid on







     
    
    
