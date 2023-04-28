clc; clear; close all

c11 = 0;
c00 = 0;
c10 = 2;
c01 = 1;

pi0 = 0.5;
pi1 = 0.5;
X = 2;
sigma = sqrt(4);

N = 1; %nb de mesures

n = 10000;


Z = zeros(N, n);
vraies_detection = zeros(N, n);

choix_jeu = menu("Choix du jeu de données", "Répartition aléatoire", "Répartition déterministe");

if choix_jeu == 1
    for i = 1:N
        nb_h0 = 0;
        nb_h1 = 0;
        for j = 1:n
            bruit = sigma * randn(1, 1);
            if (randn(1, 1) > 0 && nb_h0 < pi0*n )%choix aleatoire de l'hypothese
                Z(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            elseif (nb_h1 < pi1*n)
                Z(i, j) = bruit + X;
                nb_h1 = nb_h1 + 1;
                vraies_detection(j) = 1;
            else
                Z(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            end
        end
    end
else
    for i = 1:N
        nb_h0 = 0;
        nb_h1 = 0;
        for j = 1:n
            bruit = sigma * randn(1, 1);
            if (nb_h0 < n/2)%choix aleatoire de l'hypothese
                Z(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            elseif (nb_h1 < pi1*n)
                Z(i, j) = bruit + X(i);
                nb_h1 = nb_h1 + 1;
                vraies_detection(j) = 1;
            end
        end
    end
end

save("data", "Z", "n", "c00", "c01", "c10", "c11", "sigma", "X", "pi0", "pi1", "vraies_detection", "N");

figure(1)
plot(Z(1,:))
grid()
title("Répartition des données")