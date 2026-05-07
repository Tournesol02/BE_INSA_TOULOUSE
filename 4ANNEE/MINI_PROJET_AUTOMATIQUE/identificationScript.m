%% Initialisation
clc
clear all
close all
%load expe1.mat simout
%load expe2_multisine.mat
%load expe3_multisine.mat
%load expe4_multisine.mat
load expe5_multisine.mat

%% Variables
N = 2^9; % Nombre
Te = 5e-2; % Temps d'échantillonage
fdeb = 0; % fréquence de début d'échantillonage
ffin = 4; % fréquence de fin d'éch.
BP = [fdeb, ffin]; % Bande passante
REV = false; % Si on passe de fdeb à ffin ou l'inverse (booléen)
SHOW = true; % Montre les résultats (graphes) du MLBS (booléen)

%[u, t] = insapack.mlbs(N, Te, BP, REV, SHOW);
[u, t] = insapack.multisine(N, Te, BP, false,'all',REV, SHOW);
u = 4*u;
% var.time = t.';
% var.signals.values=u.';
% var.signals.dimensions= 1 ;
%%
% Muligens litt for hoey frekvens, den maa kanskje skrus litt ned. Deretter
% finner vi modell ved bruk av div. verktoey.

% figure()
% hold on;
% plot(simout.time, simout.signals.values)
% plot(var.time, var.signals.values)

%%% N4SID
tn      = simout.time; % On enlève les dix primières valeurs
un      = var.signals.values; % On enlève les dix primières valeurs
yn      = simout.signals.values; % On enlève les dix primières valeurs
data    = iddata(yn,un,Te);
ordre   = 2;
sys     = n4sid(data, ordre, 'Ts',Te);
y       = lsim(ss(sys), un, tn);
H       = tf(sys);
G       = (H)/(1-H);
numG = G.Numerator;
denumG = G.Denominator;


%%% Freq
[f0,U0,Y0,G0,sigU2,sigY2,sigUY2,sigG2,b,rho]    = insapack.non_param_freq(un,yn,Te);
w0                                              = 2*pi*f0;

%%% Identification via N4SID in frequency-domain
data            = iddata(Y0,U0,Te,'Frequency',w0); %Ici peut-etre changer avec f0, comme ça il marche très mieux avec Loewner
Hn4sid_fd       = n4sid(data,ordre);
y2              = lsim(ss(Hn4sid_fd), un, tn);

%% Fonction de transfert
[num, den] = ss2tf(Hn4sid_fd.a, Hn4sid_fd.b, Hn4sid_fd.c, Hn4sid_fd.d)
system_tf = tf(num, den)


%%% Identification via Loewner
wid             = 2*pi*f0(f0<BP(end)/2);
wRange          = 1:floor(length(wid)/2)*2;
wid             = 2*pi*f0(wRange);
G0id            = G0(wRange);
[la,mu,W,V,R,L] = insapack.data2loewner(wid,G0);
opt.target      = ordre;
[hr,info]       = insapack.loewner_tng(la,mu,W,V,R,L,opt);
Hloe            = dss(info.Ar,info.Br,info.Cr,info.Dr,info.Er);
Hloe            = stabsep(Hloe);
y3              = lsim(ss(Hloe), un, tn);

%%% Plots
figure()
subplot(211)
plot(tn, un)
legend("Signal d'entrée (d'exitation)")
subplot(212)
hold on;
plot(tn, yn)
plot(tn, y)
plot(tn, y2)
plot(tn, y3)
legend('Signal de sortie observée','Sortie simuléee avec n4sid en temporel', 'Sortie simuléee avec n4sid en fréquenciel', 'Sortie simuléee avec Loewner')


var2 = var;
var2.signals.values(1:100) = 1;
var2.signals.values(101:200) = 3.3;
var2.signals.values(201:300) = 0;
var2.signals.values(301:400) = 1.5;
var2.signals.values(401:503) = 2.7;

%% Comment est-ce qu'on a trouvé la valeur du P dans notre système rail
Gc=d2c(G)
consigne = 1;
P = [ 0.8 1 1.5 2 3 5];
figure()
hold on;

for i = 1:1:length(P)
    Gcf = P(i)*Gc/(1+P(i)*Gc);
    step(Gcf,10)
end;
yline(0);
title("Reponse d'échelon du système rail");
legend("Kp=0.8","Kp=1","Kp=1.5","Kp=2","Kp=3","Kp=5");


%% Bille
g = 9.81

a = 10
b = 0.2 % V/cm sur le rail, +-10V pour le rail de +-50cm
Kb = b*g
T = 0.22
Kp = 1

Gb = tf([Kb], [1 0 0])

figure()
bode(Gb)
%margin(Gb)
title('Bode système rail')

Cb = tf([a*T 1], [T, 1])*Kp; %Correcteur d'avance de phase

figure()
bode(Cb)
title('Bode correcteur avance de phase')

%% marge de phase
figure()
title("Marge de phase et gain pour le système complet");
Hc = d2c(H);
Lbo = Cb*Hc*Gb;
Lbf = Lbo / (1+Lbo);
Margin = allmargin(Lbf);
SysO=bodeplot(Lbo);
margin(Lbo);
title({'Marge de phase et gain pour le système complet BO','Gm= Inf Pm =-40.6 deg (at 5.97rad/s)'});

%% Reponse du système entier
var3.time = [0:Te:1000*Te]';
var3.signals.values = zeros(1, 1001)'
var3.signals.dimensions= 1 ;
var3.signals.values(1:150) = 0;
var3.signals.values(151:350) = -20;
var3.signals.values(351:550) = 10;
var3.signals.values(551:700) = -30;
var3.signals.values(701:1000) = 0;

Tf = var3.time(end);


load simout4_vals.mat;
load simout3_vals.mat;
pos_sortie2 = simout4.signals.values;
pos_entre = var3.signals.values(1:928,:);
commande4 = simout3.signals.values + 5; % Ajout de offset, car à plat la tension vaut +5V
Tf3 = var3.time(1:length(simout3.signals.values),:);
figure()
hold on;
Tf2 = var3.time(1:length(simout4.signals.values),:);

plot(Tf2,pos_entre);
plot(Tf2,pos_sortie2);
plot(Tf3,commande4);
legend('Position d entrée','Position de sortie','Commande en tension')
title("Reponse du système entièr");





