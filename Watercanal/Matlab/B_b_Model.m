load ModelForControl3.mat


N = length(E1);
L = 10;
a = 0;b = L;
h = (b-a)/(N-1);
zeta = a:h:b;
Internal_Damping = 0.2;

E = blkdiag(E1,E2);
Ei = E^(-1);
J = [D*0,       D;
     -D',       -Internal_Damping*E2];
B = [-phi_b,phi_a;phi_b*0,phi_a*0];

J = [D*0,D;
     -D',-Internal_Damping*E2];
B = [phi_a,-phi_b;
     phi_b*0,-phi_b*0];

rho = 1;
g = 9.8;
bb = 1;

Psi_q = 1/(2*rho)*Psi;
Psi_a = 1/(rho)*Psi;

Q = E1*rho*g/bb;


C_q = [eye(N),zeros(N)];
C_a = [zeros(N),eye(N)];



Ei  = full(Ei);
Psi = full(Psi);
J   = full(J);
B   = full(B);
E1  = full(E1);

h0 = 4;
q0      = h0*ones(N,1);
alpha0  = 0*ones(N,1);
X0  = [q0;alpha0];


%Transformation to obtain h and u from q and alpha
T = blkdiag(eye(N)/bb,eye(N)/rho);Ti = inv(T);
%TrueOutput
Ctrue = zeros(2,2*N);
Ctrue(1,1) = 1;
Ctrue(2,N) = 1;


syms X1 X2 X3 X4 X5 X6 real
syms U1 U2 real
q_d     = [X1;X2;X3];
alpha_d = [X4;X5;X6];
X = [q_d;alpha_d];
U = [U1;U2];

ed = Ei*[1/(2*rho)*Psi*kron(alpha_d,alpha_d)+E1*rho*g/bb*q_d;
         1/(rho)*Psi*kron(q_d,alpha_d)];


dXdt = Ei*(J*ed+B*U);

f = dXdt;
h = [X1;X3];

Pt_eq = solve(f == 0, X);
f_jacobian = jacobian(f,X);
%f_lin = subs(f_jacobian, [X1 X2 X3 X4 X5 X6], [1 2 3 4 5 6].');
