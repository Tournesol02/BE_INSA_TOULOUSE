function [maxreel,ecartmax]=sensible(A,pmax)
%
%  function [maxreel,ecartmax]=sensible(a,pmax)
%  Cette fonction calcule et represente dans le plan de Laplace les valeurs
%  propres de la matrice a dans le cas de perturbations non structurees
%  aleatoires, chaque element de la matrice de perturbation etant inferieur
%  ou egal a pmax. (Attention � ne pas prendre pmax trop grand aux risques
%  que la fonction confonde les modes de a)
%
%  Arguments de sortie :
%  maxreel: partie reelle maximale sur l'ensemble des valeurs propres du systeme
%  perturbe
%  ecartmax: ecart maximum relatif des valeurs propres par rapport a celles
%  du systeme nominal.
%
%  Bernard Pradin,  octobre 2005
%  revision Jeremy Leduc, Novembre 2012

% CARACTERISTIQUES DU SYSTEME NOMINAL

[n,n]=size(A);      % dimensions de a
lambda=sort(eig(A)); % Pour placer les modes dans l'ordre de leur module croissant

% GENERATION DE MATRICES PERTURBEES ET CALCUL DE VALEURS PROPRES

lw = 4;
ms = 20;
rand('seed',0);

 e=pmax*rand(n);
 p=-1;
 ensemble_lambda_perturbe=[];
 ecart_lambda=[];
 while p<= 1
    a_perturbe=A+p*e;
    lambda_perturbe=sort(eig(a_perturbe)); %Pour placer les modes dans l'ordre de leur module croissant
    ensemble_lambda_perturbe=[ensemble_lambda_perturbe;lambda_perturbe];
    ecart_lambda=[ecart_lambda;abs(lambda-lambda_perturbe) ./ abs(lambda) ];
    p=p+0.1;
 end      % while

hold off

% REPRESENTATION DES VALEURS PROPRES NOMINALES

% Definition du cadre graphique et trace premiere val. prop.

x1=max(0,max(real(ensemble_lambda_perturbe))+0.1);
axis([1.1*min(real(ensemble_lambda_perturbe))-1,x1,1.1*min(imag(ensemble_lambda_perturbe))-1,1.1*max(imag(ensemble_lambda_perturbe))+1])
plot(real(lambda(1)),imag(lambda(1)),'rh','LineWidth',lw,'MarkerSize',ms)
grid

hold on

% Trace autres valeurs propres

for i=2:n
     plot(real(lambda(i)),imag(lambda(i)),'rh','LineWidth',lw,'MarkerSize',ms)
end

% REPRESENTATION DES VALEURS PROPRES DE LA MATRICE PERTURBEE

nm=max(size(ensemble_lambda_perturbe));

for i=1:nm
    plot( real(ensemble_lambda_perturbe(i)),imag(ensemble_lambda_perturbe(i)),'.','LineWidth',lw,'MarkerSize',ms)
end

set(gca,'FontSize',20)
hold off

% CALCUL DES ARGUMENTS DE SORTIE

maxreel=max(real(ensemble_lambda_perturbe));
ecartmax=max(ecart_lambda);

