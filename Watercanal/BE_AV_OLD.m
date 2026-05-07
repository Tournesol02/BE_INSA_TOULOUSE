m = 100000;
g = 9.8;
c = 3;
Kc = 1500;
I = 64*m;
Gamma = 192*m;
d1 = 1;
d2 = 25;
h = 0.5;
t = 0:0.2:100;

X1_eq = 1000/3.6;
U1_eq = c*X1_eq^2; % On met climb angle X2 à 0
X2_eq = 0;
X3_eq = (m*g/X1_eq) / (Kc*X1_eq + U1_eq/X1_eq);
X4_eq = 0;
U2_eq = (Kc*d1*X1_eq^2*X3_eq - U1_eq*h) / d2;

X_eq = [X1_eq; X2_eq; X3_eq; X4_eq];
U_eq = [U1_eq ; U2_eq];

[t,X] = ode45(MyODE(t,X_eq,U_eq), t, X_eq);



h1 = X1_eq;
h2 = X1_eq*sin(X2_eq);

function dXdt = MyODE(t,X, U)
X1 = X(1);
X2 = X(2);
X3 = X(3);
X4 = X(4);
U1 = U(1);
U2 = U(2);

m = 100000;
g = 9.8;
c = 3;
Kc = 1500;
I = 64*m;
Gamma = 192*m;
d1 = 1;
d2 = 25;
h = 0.5;

d_X1 = (1/m)*(U1 - m*g*sin(X2) - c*X1^2);
d_X2 = (1/m)*(Kc * X1 *(X3-X2) + U1*(X3-X2)/X1 - m*g*cos(X2)/X1);
d_X3 = X4;
d_X4 = (1/I) * (-Gamma*X4 - Kc*d1*X1^2 * (X3-X2) + U2*d2 + U1*h);

d_X = [d_X1; d_X2; d_X3; d_X4];
dXdt = d_X;
end
