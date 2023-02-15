clc; clear; close all

c11 = 0;
c00 = 0;
c10 = 2;
c01 = 1;
pi0 = 0.5;
pi1 = 0.5;
v = 5;
sigma = 3;

m = 1; %nb de mesures

N = 10000;


X = zeros(m, N);

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



fonction = exp(-(1/2)*((X-v)/sigma).^2) .* exp((1/2)*(X/sigma).^2);
lambda_0 = (pi0*(c10-c00))/(pi1*(c01-c11)); %seuil detecteur de bayes
seuil = (2*sigma^2*log(lambda_0) + v^2)/(2*v); %seuil detecteur equivalent

d1 = 1;
d0 = 0;

%Detecteur de Bayes
delta1 = zeros(1, N);
compteur_Bayes = 0;
for i = 1:m
    for j = 1:N
        if (fonction(j) > lambda_0)
            delta1(j) = d1;
            compteur_Bayes = compteur_Bayes + 1;
        else
            delta1(j) = d0;
        end
    end
end

%Detecteur equivalent
delta2 = zeros(1, N);
compteur_equivalent = 0;
for i = 1:m
    for j = 1:N
        if (X(i,j) > seuil)
            delta2(j) = d1;
            compteur_equivalent = compteur_equivalent + 1;
        else
            delta2(j) = d0;
        end
    end
end

for i = 1:m
    nb_d1_detecte = 0;
    nb_fa = 0;
    for j = 1:N
        if (delta1(j) == 1 && vraies_detection(j) == 1)
            nb_d1_detecte = nb_d1_detecte + 1;
        elseif (delta1(j) == 1 && vraies_detection(j) == 0)
            nb_fa = nb_fa + 1;
        end
    end
end
Pd = nb_d1_detecte / (pi1*N)
Pfa = nb_fa/(pi0*N)

figure(1)
plot(X(1,:))
grid()
title("Mesures X")