

close all
clear all
clc

%Parameters
m       = 100000;
g       = 9.8;
Kc      = 1500;
c       = 3;
I       = 64*m;
Gamma   = 192*m;
d1_bar  = 1;
d2_bar  = 25;
h_bar   = 0.5;

%Save the parameters in a column vector for ode45 simulation
Parameters =[m;g;Kc;c;I;Gamma;d1_bar;d2_bar;h_bar];


%We set the equilibrium point
eta_eq      = 0;
omega_eq    = 0;
X1_eq       = 1000/3.6;
% T_eq        = 90000;

X2_eq = eta_eq;
X4_eq = omega_eq;
% U1_eq = T_eq;

% X1_eq = sqrt(U1_eq/c);
U1_eq = c*X1_eq^2;
X3_eq = m*g/(Kc*X1_eq^2+U1_eq);
U2_eq = (Kc*d1_bar*X1_eq^2*X3_eq-U1_eq*h_bar)/d2_bar;

X_eq = [X1_eq;X2_eq;X3_eq;X4_eq];
U_eq = [U1_eq;U2_eq];
Y_eq = [X1_eq;X1_eq*sin(X2_eq)];

%Verify Equilibrium
f1_eq = 1/m*(U1_eq-m*g*sin(X2_eq)-c*X1_eq^2);
f2_eq = 1/m*(Kc*X1_eq*(X3_eq-X2_eq) + U1_eq*(X3_eq-X2_eq)/X1_eq - m*g*cos(X2_eq)/X1_eq);
f3_eq = X4_eq;
f4_eq = 1/I*(-Gamma*X4_eq - Kc*d1_bar*X1_eq^2*(X3_eq-X2_eq) + U2_eq*d2_bar + U1_eq*h_bar) ;


%Linearization
syms X1 X2 X3 X4 real
syms U1 U2 real

f1 = 1/m*(U1-m*g*sin(X2)-c*X1^2);
f2 = 1/m*(Kc*X1*(X3-X2) + U1*(X3-X2)/X1 - m*g*cos(X2)/X1);
f3 = X4;
f4 = 1/I*(-Gamma*X4 - Kc*d1_bar*X1^2*(X3-X2) + U2*d2_bar + U1*h_bar) ;

h1 = X1;
h2 = X1*sin(X2);

f = [f1;f2;f3;f4];
h = [h1;h2];
X = [X1;X2;X3;X4];
U = [U1;U2];

%Linearization arround the equilibirum
%Linearize the nonlinear equations around the equilibrium f(X_eq,U_eq) = 0.
A = jacobian(f,X);
B = jacobian(f,U);
C = jacobian(h,X);
D = jacobian(h,U);


X1 = X_eq(1);
X2 = X_eq(2);
X3 = X_eq(3);
X4 = X_eq(4);

U1 = U_eq(1);
U2 = U_eq(2);

%% Question 2
fprintf('Question 2')

X_eq
U_eq
Y_eq


A = eval(A)
B = eval(B)
C = eval(C)
D = eval(D)


%% Question 3
SYS = ss(A,B,C,D);

eig(SYS)

%% Question 4
fprintf("Le système est observable et contrôlable car le rang des matrices Qo et Qc est égal à l'ordre du SYStème (=4).")
rank(obsv(SYS))
rank(ctrb(SYS))

fprintf('For the stability analysis, we use eig(SYS). \n\n')
eig(SYS)


%% Q5
fprintf('The Matlab function damp(SYS) shows the eigenvalues, the damping coefficients and the natural frequencies. \n\n')
damp(SYS)
% pause

lw = 4;
ms = 20;
fs = 16;
fprintf('\n The Matlab function pzmap(SYS) shows in a figure the poles and zeros of SYS. \n\n')
[P,Z] = pzmap(SYS);
figure
pzmap(SYS)
hold on
plot(real(P),imag(P),'x','LineWidth',lw,'MarkerSize',ms)
plot(real(Z),imag(Z),'o','LineWidth',lw,'MarkerSize',ms)
grid
title('Poles et zeros du SYSteme dans le plan de Laplace','Interpreter','latex','FontSize',fs)
% xlabel('Real axis','Interpreter','latex','FontSize',fs)
% ylabel('Imaginary axis','Interpreter','latex','FontSize',fs)
% set(gca,'FontSize',fs*0.8)

TFINAL = 200;
fprintf('\n The Matlab function impulse(SYS) shows the impulse response. \n\n')
figure
impulse(SYS,TFINAL)
% pause


fprintf('\n The Matlab function step(SYS) shows the step response. \n\n')
figure
step(SYS,TFINAL)



%% Simulation
%Reference signals
r1 = @(t) 1*U_eq(1)*heaviside(t-5);
r2 = @(t) 1*U_eq(2)*(heaviside(t-20) - 2*heaviside(t-40));
% r2 = @(t) 1*heaviside(t-10);


%Simulation of the nonlinear SYStem
X_0     = X_eq;   %Initial condition at the equilibrium
dt      = 0.2;                              %Sample time
t       = 0:dt:100;
[t,X]   = ode45(@(t,X) MyODE(t,X,Parameters,U_eq,r1,r2),t,X_0);
t       = t';
X       = X';
V       = X(1,:);
eta     = X(2,:);
theta   = X(3,:);
omega   = X(4,:);


U1 = r1(t) + U_eq(1);
U2 = r2(t) + U_eq(2);

Y1 = V;
Y2 = V.*sin(eta);


%Simulation of the linearized SYStem

%x is the linearized variable x = X - X_eq
x_0     = X_0 - X_eq;
u       = [U1'- U_eq(1),U2' - U_eq(2)];
[y,t,x] = lsim(SYS,u,t,x_0); 
x       = x';
y       = y';
t       = t';


V_lin      = x(1,:) + X_eq(1);
eta_lin    = x(2,:) + X_eq(2);
theta_lin  = x(3,:) + X_eq(3);
omega_lin  = x(4,:) + X_eq(4);


y1 = y(1,:);
y2 = y(2,:);

%%
PlotParameters
figure
set(gcf,'units','points','position',[x0screen,y0screen,WidthScreen,HeightScreen])
subplot(4,1,1)
title({'States'},'Interpreter','latex','FontSize',fs*4)
hold on
yline(X_eq(1),'--','LineWidth',lw/2);
plot(t,V,'LineWidth',lw)
plot(t,V_lin,'LineWidth',2)
legend({'$X_1^\star$','$X_1(t)$','$x_1(t) + X_1^\star$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)

subplot(4,1,2)
hold on
yline(X_eq(2),'--');
plot(t,eta,'LineWidth',2)
plot(t,eta_lin,'LineWidth',2)
legend({'$X_2^\star$','$X_2(t)$','$x_2(t) + X_2^\star$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)

subplot(4,1,3)
hold on
yline(X_eq(3),'--');
plot(t,theta,'LineWidth',2)
plot(t,theta_lin,'LineWidth',2)
legend({'$X_3^\star$','$X_3(t)$','$x_3(t) + X_3^\star$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)

subplot(4,1,4)
hold on
yline(X_eq(4),'--');
plot(t,omega,'LineWidth',2)
plot(t,omega_lin,'LineWidth',2)
legend({'$X_4^\star$','$X_4(t)$','$x_4(t) + X_4^\star$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)
xlabel({'time $t$ $(s)$'},'Interpreter','latex','FontSize',fs)



figure
set(gcf,'units','points','position',[x0screen,y0screen,WidthScreen,HeightScreen])
subplot(2,2,1)
hold on
yline(U_eq(1),'--','LineWidth',lw/2);
plot(t,U1,'LineWidth',lw)
legend({'$U_1^\star$','$U_1(t)$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)

subplot(2,2,2)
hold on
yline(Y_eq(1),'--','LineWidth',lw/2);
plot(t,Y1,'LineWidth',lw)
plot(t,y1+Y_eq(1),'LineWidth',lw)
legend({'$Y_1^\star$','$Y_1(t)$','$y_1(t) + Y_1^\star$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)

subplot(2,2,3)
hold on
yline(U_eq(2),'--');
plot(t,U2,'LineWidth',2)
legend({'$U_2^\star$','$U_2(t)$'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)
xlabel({'time $t$ $(s)$'},'Interpreter','latex','FontSize',fs)


subplot(2,2,4)
hold on
yline(Y_eq(2),'--','LineWidth',lw/2);
plot(t,Y2,'LineWidth',lw)
plot(t,y2+Y_eq(2),'LineWidth',lw)
legend({'$Y_2^\star$','$Y_2(t)$','$y_2(t)+Y_2^\star $'},'Interpreter','latex','FontSize',fs,'Location','Northwest')
set(gca,'FontSize',fst)



pause(2)

Animation(V,omega,t)

%%

% Question 8
% U est de la taille 2*1, X est de la taille 4*1, donc K est forcément est
% de la taille 2*4. Pour H, R est 2*1, donc il est de la taille 2*2.

%% Question 9
% Le temps de réponse est à 3*tau, et on veut un temps de réponse à 30 sec.
% Donc /lambda1 = -3/30 = -0.1. \lambda2 = -3/20 car Tr = 20s.
% \lambda3 = -5 et \lambda4 = - 10.
lambda_1 = -3/30;
lambda_2 = -3/20;

%% Question 10
% La fonction place prend en arguments les matrices A et B, et les pôles
% souhaités. Elle renvoie la matrice K dans le cas d'un retour d'état X =
% -K*U.
% La fonction sensible prend en argument la matrice A en boucle fermée et
% donne à quel point les valeurs popres voulues sont déplacées dans le plan
% des racines

%% Question 11
% Puisque le système est commandable, on peut assigner lespôles souhaités.
poles = [lambda_1, lambda_2, -10, -30];
K_1 = place(A,B, poles)

%% Question 12
rank(ctrb(A,B(:,1)))
K_2 = place(A,B(:,1), poles)

%% Question 13
rank(ctrb(A,B(:,2)))
K_3 = place(A,B(:,2), poles)

%% Question 14
[~,~] = sensible(A-B*K_1, 0.2);
figure
[~,~] = sensible(A-B(:,1)*K_2, 0.2);
figure
[~,~] = sensible(A-B(:,2)*K_3, 0.2);

%% Question 15
H = inv(C*inv(-A+B*K)*B);
sysBF = ss((A-B*K), B*H, C, D)
step(sysBF)
%%

function dXdt = MyODE(t,X,Parameters,U_eq,r1,r2)
X1   = X(1);
X2   = X(2);
X3   = X(3);
X4   = X(4);

m       = Parameters(1);
g       = Parameters(2);
Kc      = Parameters(3);
c       = Parameters(4);
I       = Parameters(5);
Gamma   = Parameters(6);
d1_bar  = Parameters(7);
d2_bar  = Parameters(8);
h_bar   = Parameters(9);

U1 = r1(t) + U_eq(1);
U2 = r2(t) + U_eq(2);


f1 = 1/m*(U1-m*g*sin(X2)-c*X1^2);
f2 = 1/m*(Kc*X1*(X3-X2) + U1*(X3-X2)/X1 - m*g*cos(X2)/X1);
f3 = X4;
f4 = 1/I*(-Gamma*X4 - Kc*d1_bar*X1^2*(X3-X2) + U2*d2_bar + U1*h_bar) ;

f = [f1;f2;f3;f4];

dXdt = f;
end 


%Q1)
