clc; clear; close all

c11 = 0;
c00 = 0;
c10 = 2;
c01 = 1;

pi0 = 0.5;
pi1 = 0.5;
v = 2;
sigma = sqrt(4);

m = 1; %nb de mesures

N = 10000;


X = zeros(m, N);
vraies_detection = zeros(N, 1);

choix_jeu = menu("Choix du jeu de données", "Répartition aléatoire", "Répartition déterministe");

if choix_jeu == 1
    for i = 1:m
        nb_h0 = 0;
        nb_h1 = 0;
        for j = 1:N
            bruit = sigma * randn(1, 1);
            if (randn(1, 1) > 0 && nb_h0 < pi0*N )%choix aleatoire de l'hypothese
                X(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            elseif (nb_h1 < pi1*N)
                X(i, j) = bruit + v;
                nb_h1 = nb_h1 + 1;
                vraies_detection(j) = 1;
            else
                X(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            end
        end
    end
else
    for i = 1:m
        nb_h0 = 0;
        nb_h1 = 0;
        for j = 1:N
            bruit = sigma * randn(1, 1);
            if (nb_h0 < N/2)%choix aleatoire de l'hypothese
                X(i, j) = bruit;
                nb_h0 = nb_h0 + 1;
                vraies_detection(j) = 0;
            elseif (nb_h1 < pi1*N)
                X(i, j) = bruit + v;
                nb_h1 = nb_h1 + 1;
                vraies_detection(j) = 1;
            end
        end
    end
end

save("data", "X", "c00", "c01", "c10", "c11", "sigma", "v", "pi0", "pi1", "vraies_detection");

figure(1)
plot(X)
grid()
title("Répartition des données")