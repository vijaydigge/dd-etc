clear all;
n = 2;
T = 20;

%% Data collection
% A = [1.01 0.0995; 0.01 0.9900];
% 
% B = [0.01; -0.184];

A = [0.91 0.0995; 0.02 0.9900];
B = [0.01; -0.184];

% rng(1,'philox');    /      %random input 'u1' is same all the time
% u1=1*(rand(T,1));

% for k = 1:T
% x(:,k+1) = A*x(:,k) + B*[u1(k)];
% end

pol = [];
si = [];
jj = [];
for j = 0.001: 0.005: 0.5
% varr = 0.00000010;
jj(end+1) = j;
rng(1,'philox');          %random input 'u1' is same all the time
u1=1*(rand(T,1));
x(:,1) = [0.5; -0.4];
for k = 1:T
x(:,k+1) = A*x(:,k) + B*[u1(k)];
end
x = x+ j*(randn(2,T+1)-randn(2,T+1));
u = u1';
Z = [u; x(:,1:T)];

X1 = x(:,2:T+1);
X0 = x(:,1:T);
% 
% figure(1);
% plot(u);
% hold on;
% plot(X0(1,:));
% hold on;
% plot(X0(2,:));
%%  Stability Analysis
%  load('data_et.mat')  %system 1

%  load('et-data_batch_react.mat') % system 2

lam = 0.85;   %signifieig(A+B*K)es how quickly states have to converge
cvx_begin sdp
    variables Q1(T,n) ;
%     minimize trace(X0*Q1)
%     maximize trace(X0*Q1) 
     %P == inv(X0*Q);
    [lam*X0*Q1 Q1'*X1' ;  X1*Q1 X0*Q1] >= 0*eye(Noise order2*n);
    X0*Q1 +  (X0*Q1)' >= 0.00001*eye(n,n);
     cvx_end
 K = u*Q1*inv(X0*Q1);
 
%%
% rng(1,'philox');          %random input 'u1' is same all the time
% u1=(rand(T,1));
% 
% for k=1:T
% x(:,k+1) = A*x(:,k) + B*[u1(k)];
% end


%% ET 
q = 7;
 cvx_begin sdp
    variable q;
    variable sig;
    variable Q2(T,n);
    maximize sig
    subject to
        q<=7.5;
        sig >= 0.0001;
        sig <= 0.9999;
        q >= 0.001;      
       [lam*X0*Q1 0*eye(n) Q1'*X1' sig*X0*Q1; 0*eye(n) q*eye(n) Q2'*X1' 0*eye(n); X1*Q1 X1*Q2 X0*Q1 0*eye(n); sig*(X0*Q1) 0*eye(n) 0*eye(n) q*eye(n)] >= 0.00*eye(4*n);   
         u*Q2 - q*K == 0 ;Noise order
        X0*Q2 == 0*eye(n);
  cvx_end
  

flag = 0;
si(end+1) = sig;
pol(:,end+1) = max(abs(eig(A+B*K)));
% if(pol(1,:)>=1)
%     flag = flag+1;
% end
end
figure(2);
ang=0:0.01:2*pi; 
xp=1*cos(ang);
yp=1*sin(ang);
xpla=sqrt(lam)*cos(ang);
ypla=sqrt(lam)*sin(ang);
hold on
plot(xp,yp,'r','LineWidth',1.5);
plot(xpla,ypla,'b','LineWidth',1.5);
stem(pol(1,:),zeros(length(pol),1),'LineWidth',1.5);
% stem(pol(2,:),zeros(length(pol),1));
hold off;

figure(1);
hold on;
plot(0.001: 0.005: 0.5,si,'k','LineWidth',1.5);
plot(0.001: 0.005: 0.5,pol(1,:),'r','LineWidth',1.5);
plot(0.001: 0.005: 0.5,sqrt(0.85)*ones(1,length(si)),'b','LineWidth',1.5);
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Noise level}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{Magnitude}','Interpreter','Latex','fontsize',22,'fontweight','bold');
h1=legend('\bf{$\sigma$}','\bf{$\lambda_{max}(X_+G)$}','\bf{$\sqrt{\lambda}$}');
set(h1,'fontsize',22,'Interpreter','Latex')
box on;
xline(0.131,'LineWidth',1.5);
% figure(1);
% ang=0:0.01:2*pi; 
% xp=1*cos(ang);
% yp=1*sin(ang);
% xpla=sqrt(lam)*cos(ang);
% ypla=sqrt(lam)*sin(ang);
% hold on
% plot(xp,yp,'r');
% plot(xpla,ypla,'b');
% stem(pol(1,:),zeros(100,1));
% stem(pol(2,:),zeros(100,1));
% hold off;

% P = inv(X0*Q1);