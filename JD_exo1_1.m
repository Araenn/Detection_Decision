clc; clear; close all

c11 = 0;
c00 = 0;
c10 = 2;
c01 = 1;

pi0 = 0.5;
pi1 = 0.5;
sigma = sqrt(4);

N = 1; %nb de mesures

n = 10000;


Z = zeros(N, n);
vraies_detection = zeros(N, n);

choix_jeu = menu("Choix du jeu de données", "X = 2", "X = 4", "X = 8", "X = 10");

if choix_jeu == 1
    X = 2;
elseif choix_jeu == 2
    X = 4;
elseif choix_jeu == 3
    X = 8;
else
    X = 10;
end
nb_h0 = 0;
nb_h1 = 0;
for i = 1:n
    bruit = sigma * randn(1, 1);
    if (nb_h0 < n/2)%choix aleatoire de l'hypothese
        Z(i) = bruit;
        nb_h0 = nb_h0 + 1;
        vraies_detection(i) = 0;
    elseif (nb_h1 < pi1*n)
        Z(i) = bruit + X;
        nb_h1 = nb_h1 + 1;
        vraies_detection(i) = 1;
    end
end

save("data", "Z", "n", "c00", "c01", "c10", "c11", "sigma", "X", "pi0", "pi1", "vraies_detection", "N");

figure(1)
plot(Z(1,:))
grid()
title("Répartition des données")