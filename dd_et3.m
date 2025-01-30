
clear all;
n = 2;

T = 20;
%% Data collection
% A = [0.1 1.2; 0.007 1.05];
% A = [1.0018 0.01; 0.36 1.0018]
% A = [1.01 2; 0.01 0.9]
%  A = [0.1 1.2;0.007 2.05];
%  B = [300 200;0.5 0.001];
%  B = [0.03; 0.05];
 
% B = [-0.01; -1.84];
%  A = [1.178 0.001; -0.051 0.661];
%  B = [-0.004;-0.467];

A = [0.9 0.0995; 0.02 0.9900];
B = [0.01; -0.184];

x(:,1) = [0.5; -0.4];
rng(1,'philox');          %random input 'u1' will be same all the time
u1=01*(rand(T,1));
% u2=01*(rand(T,1));
% u = [u1';u2'];
u = u1';
for k=1:T
%     x(:,k+1) = A*x(:,k) + B*[u1(k);u2(k)];
x(:,k+1) = A*x(:,k) + B*[u(k)];
end


% Z = [u; x(:,1:T)];

% rank(Z); %%persistently exciting data condition
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

lam = 0.85;   %signifies how quickly states have to converge
cvx_begin sdp
    variables Q1(T,n) ;
%     minimize trace(X0*Q1)
%     maximize trace(X0*Q1) 
     %P == inv(X0*Q);
    [lam*X0*Q1 Q1'*X1' ;  X1*Q1 X0*Q1] >= 0*eye(2*n);
    X0*Q1 +  (X0*Q1)' >= 0.00001*eye(n,n);
%     X0*Q1 <= 1*eye(n);
     cvx_end
 K = u*Q1*inv(X0*Q1);
 %% Eigen value check
eig(A+B*K)

%% ET 
% q = 8
 cvx_begin sdp
    variables  sig
    variable q
    variable Q2(T,n)
    maximize sig
    subject to
        sig >= 0.0001;
        sig <= 0.9999;
        q >= 0.0001;      
       [lam*X0*Q1 0*eye(n) Q1'*X1' sig*X0*Q1; 0*eye(n) q*eye(n) Q2'*X1' 0*eye(n); X1*Q1 X1*Q2 X0*Q1 0*eye(n); sig*(X0*Q1) 0*eye(n) 0*eye(n) q*eye(n)] >= 0.00*eye(4*n);   
         u*Q2 - q*K == 0 ;
        X0*Q2 == 0*eye(n);
  cvx_end
  
  
%   A = [0.1 1.2; 0.007 1.05];  %system 1
% B = [300 200; 0.5 0.001];

% A = [1.178 0.001 0.511 -0.403; -0.051 0.661 -0.011 0.061; 0.076 0.335 0.560 0.382; 0 0.335 0.089 0.849];
% B = [0.004 -0.087; 0.467 0.01; 0.213 -0.235; 0.213 -0.016];   %system 2

%% Plots
sigd = cvx_optval; %0.283465
% sigd = 0.07
v = [];
flag=[];
e=[];
xhat=x;
kmax = 50;
for i = 1:2*kmax
    e(:,i) = xhat(:,i) - x(:,i);
    norme(i) = norm(e(:,i));
      normx(1) = norm(x);
if norm(e(:,i))< sigd*norm(x(:,i))
             x(:,i+1) = A*x(:,i)+B*K*xhat(:,i);
             normx(i+1)=norm(x(:,i+1));
             if norm(xhat(:,i) - x(:,i+1)) < sigd*norm(x(:,i+1))
             xhat(:,i+1)=xhat(:,i);
             flag(i)=0;

             else 
             xhat(:,i+1)=x(:,i+1); 
             flag(i)=1;

             end
             
else
        xhat(:,i)=x(:,i);
        x(:,i+1)=A*x(:,i)+B*K*xhat(:,i);
        normx(i+1)=norm(x(:,i+1));
       %  xhat(:,i+1)=x(:,i+1);
        xhat(:,i+1)=xhat(:,i);
         flag(i)=1;
end
    V1(i)=x(:,i)'*inv(X0*Q1)*x(:,i);
    V2(i) = x(:,i+1)'*inv(X0*Q1)*x(:,i+1);
    v(i) = V2(i)/V1(i);
    %rand1(i+1)=mag*(rand(1)); 
end


indices=find(flag==1);
flag_counter=[];
for j=1:length(indices)-1
  flag_counter(j)=indices(j+1)-indices(j);
end

disp("Max IET");
max(flag_counter)
disp("Min IET");
min(flag_counter)
disp("Avg IET");
count = 0;
for i=1:length(flag_counter)
count = count+flag_counter(i);
end
avg = ceil(count/length(flag_counter))


flag_count=[];
flag_current=0;
for l=1:length(flag)
    if flag(l)== 0
        flag_current=flag_current+1;
        flag_count(l)=NaN;
    else
        flag_current=flag_current+1;
        flag_count(l)=flag_current;
        flag_current=0;
    end
end
    
figure(2);
stem(1:length(flag),flag_count, 'k','LineWidth',1.5);
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{Inter event times}','Interpreter','Latex','fontsize',22,'fontweight','bold');


figure(3);
% Box on;
hold on;
stem(1:length(norme),norme,'r','LineWidth',1.5);
plot(1:length(normx),sigd*normx,'k','LineWidth',1.5);
hold off;
h1=legend('\bf{$\|e(k)\|$}','\bf{$\sigma\|x(k)\|$}');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$\|e(k)\|$} and \bf{$\sigma\|x(k)\|$}','Interpreter','Latex','fontsize',22,'fontweight','bold');
box on
% title('Switching Law');

figure(4);
hold on;
plot(1:length(x(1,:)),x(1,:),'b-.','LineWidth',1.5);
plot(1:length(x(2,:)),x(2,:),'r-.','LineWidth',1.5);
hold off;
h1=legend('\bf{$x_1(k)$}','\bf{$x_2(k)$}');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$x_1(k)$} and \bf{$x_2(k)$}','Interpreter','Latex','fontsize',22,'fontweight','bold');
box on;

figure(5);
plot(1:length(v(1,:)), (ones(length(v),1)*lam)', 'k','LineWidth',1.5);
hold on;
stairs(1:length(v(1,:)), v,'r','LineWidth',1.5);
h1=legend('\bf{$\lambda$}','\bf{$\frac{V(k+1)}{V(k)}$}');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$\frac{V(k+1)}{V(k)}$} and \bf{$\lambda$}','Interpreter','Latex','fontsize',22,'fontweight','bold');


%% Periodic implementation
xc = [];   
normxc = [];
xc = x(:,1);
normxc = norm(xc(:,1));
normx(1) = normxc(1);
for k=1:2*kmax
xc(:,k+1) = (A+B*K)*xc(:,k);
normxc(k+1) = norm(xc(:,k+1));

 Vv1(k)=xc(:,k)'*inv(X0*Q1)*xc(:,k);
    Vv2(k) = xc(:,k+1)'*inv(X0*Q1)*xc(:,k+1);
    vv(k) = Vv2(k)/Vv1(k);
    
end
figure(6);
hold on;
plot(1:length(x(1,:)),xc(1,:),'k','LineWidth',1.5);
plot(1:length(x(1,:)),x(1,:),'r','LineWidth',1.5);
plot(1:length(x(2,:)),xc(2,:),'k','LineWidth',1.5);
plot(1:length(x(2,:)),x(2,:),'r','LineWidth',1.5);
hold off;
h1=legend('Periodic control','Event-triggered control');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$x_1(k)$} \bf{and} \bf{$x_2(k)$}','Interpreter','Latex','fontsize',22,'fontweight','bold');
box on;


%% norm of states ET vs PC
figure(7);
hold on;
plot(1:length(xc(1,:)), normxc,'k','LineWidth',1.5);
plot(1:length(x(1,:)),normx,'r','LineWidth',1.5);
hold off;
h1=legend('Periodic control','Event-triggered control');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$\|x(k)\|$}','Interpreter','Latex','fontsize',22,'fontweight','bold');
box on


figure(8);
plot(1:length(vv(1,:)), (ones(length(vv),1)*lam)', 'k','LineWidth',1.5);
hold on;
stairs(1:length(vv(1,:)), vv,'b','LineWidth',1.5);
h1=legend('\bf{$\lambda$}','\bf{$\frac{V(k+1)}{V(k)}$}');
set(h1,'fontsize',22,'Interpreter','Latex')
set(gca, 'fontsize',22, 'fontweight','bold')
xlabel('\bf{Sample}','Interpreter','Latex','fontsize',22,'fontweight','bold'); 
ylabel('\bf{$\frac{V(k+1)}{V(k)}$} and \bf{$\lambda$}','Interpreter','Latex','fontsize',22,'fontweight','bold');

