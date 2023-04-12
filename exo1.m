clc; clear; close all

load("data.mat")

fonction = exp(-(1/2)*((X-v)/sigma).^2) .* exp((1/2)*(X/sigma).^2); %detecteur de bayes
lambda_0 = (pi0*(c10-c00))/(pi1*(c01-c11)); %seuil detecteur de bayes
seuil = (2*sigma^2*log(lambda_0) + v^2)/(2*v); %seuil detecteur equivalent

d1 = 1;
d0 = 0;

%Detecteur de Bayes
delta1 = zeros(m, N);
compteur_Bayes = 0;
for i = 1:m
    for j = 1:N
        if (fonction(i, j) > lambda_0)
            delta1(i, j) = d1;
            compteur_Bayes = compteur_Bayes + 1;
        else
            delta1(i, j) = d0;
        end
    end
end

%Detecteur equivalent
delta2 = zeros(m, N);
compteur_equivalent = 0;
for i = 1:m
    for j = 1:N
        if (X(i,j) > seuil)
            delta2(i, j) = d1;
            compteur_equivalent = compteur_equivalent + 1;
        else
            delta2(i, j) = d0;
        end
    end
end

for i = 1:m
    nb_d1_detecte_bayes = 0;
    nb_fa_bayes = 0;
    nb_d1_detecte_equivalent = 0;
    nb_fa_equivalent = 0;
    
    for j = 1:N
        if (delta1(i, j) == 1 && vraies_detection(i, j) == 1)
            nb_d1_detecte_bayes = nb_d1_detecte_bayes + 1;
        elseif (delta1(i, j) == 1 && vraies_detection(i, j) == 0)
            nb_fa_bayes = nb_fa_bayes + 1;
        end
        if (delta2(i, j) == 1 && vraies_detection(i, j) == 1)
            nb_d1_detecte_equivalent = nb_d1_detecte_equivalent + 1;
        elseif (delta2(i, j) == 1 && vraies_detection(i, j) == 0)
            nb_fa_equivalent = nb_fa_equivalent + 1;
        end
    end
    
    Pfa_bayes(i) = nb_fa_bayes/(pi0*N);
    Pd_bayes(i) = nb_d1_detecte_bayes / (pi1*N);
    
    Pfa_equivalent(i) = nb_fa_equivalent/(pi0*N);
    Pd_equivalent(i) = nb_d1_detecte_equivalent / (pi1*N);
end


figure(1)
plot(X)
grid()
title("Mesures X")

figure(2)
stem(Pfa_bayes, Pd_bayes, "LineStyle", "none")
hold on
stem(Pfa_equivalent, Pd_equivalent, "LineStyle", "none")
grid()
xlabel("Pfa")
ylabel("Pd")
title("Courbes COR pour Bayes et équivalent")