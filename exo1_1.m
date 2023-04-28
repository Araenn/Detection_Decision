clc; clear; close all

load("data.mat")

fonction = exp(-(1/2)*((Z-X)/sigma).^2) .* exp((1/2)*(Z/sigma).^2); %detecteur de bayes
lambda_0 = (pi0*(c10-c00))/(pi1*(c01-c11)); %seuil detecteur de bayes
seuil = (2*sigma^2*log(lambda_0) + X.^2)/(2*X); %seuil detecteur equivalent
seuil_min_bayes = min(Z);
seuil_max_bayes = max(Z);
nb_seuils = 100;
seuils_bayes = linspace(seuil_min_bayes, seuil_max_bayes, nb_seuils);

seuil_min_eq = min(Z);
seuil_max_eq = max(Z);
seuils_eq = linspace(seuil_min_eq, seuil_max_eq, nb_seuils);



d1 = 1;
d0 = 0;

%Detecteur de Bayes
delta1 = zeros(nb_seuils, n);
compteur_Bayes = 0;
for i = 1:nb_seuils
    for j = 1:n
        if (fonction(j) > seuils_bayes(i))
            delta1(i, j) = d1;
            compteur_Bayes = compteur_Bayes + 1;
        else
            delta1(i, j) = d0;
        end
    end
end

%Detecteur equivalent
delta2 = zeros(nb_seuils, n);
compteur_equivalent = 0;
for a = 1:nb_seuils
    for b = 1:n
        if (Z(b) > seuils_eq(a))
            delta2(a, b) = d1;
            compteur_equivalent = compteur_equivalent + 1;
        else
            delta2(a, b) = d0;
        end
    end
end


for a = 1:nb_seuils
    
    nb_d1_detecte_bayes = 0;
    nb_fa_bayes = 0;
    nb_d1_detecte_equivalent = 0;
    nb_fa_equivalent = 0;
    for b = 1:n
        if (delta1(a, b) == 1 && vraies_detection(b) == 1)
            nb_d1_detecte_bayes = nb_d1_detecte_bayes + 1;
        elseif (delta1(a, b) == 1 && vraies_detection(b) == 0)
            nb_fa_bayes = nb_fa_bayes + 1;
        end
        if (delta2(a, b) == 1 && vraies_detection(b) == 1)
            nb_d1_detecte_equivalent = nb_d1_detecte_equivalent + 1;
        elseif (delta2(a, b) == 1 && vraies_detection(b) == 0)
            nb_fa_equivalent = nb_fa_equivalent + 1;
        end
    end
    
    Pfa_bayes(a) = nb_fa_bayes/(pi0*n);
    Pd_bayes(a) = nb_d1_detecte_bayes / (pi1*n);
    risque_bayes(a) = c10 * pi0 * Pfa_bayes(a) + c01 * pi1 * (1 - Pd_bayes(a) );
    
    Pfa_equivalent(a) = nb_fa_equivalent/(pi0*n);
    Pd_equivalent(a) = nb_d1_detecte_equivalent / (pi1*n);
    risque_equivalent(a) = c10 * pi0 * Pfa_equivalent(a) + c01 * pi1 * (1 - Pd_equivalent(a) );
end


figure(1)
plot(Z(1,:))
grid()
title("Mesures Z")

figure(2)
plot(Pfa_bayes, Pd_bayes)
hold on
plot(Pfa_equivalent, Pd_equivalent)
grid()
xlabel("Pfa")
ylabel("Pd")
title("Courbes COR pour Bayes et équivalent")
legend("Détecteur de Bayes", "Détecteur équivalent", 'Location', 'southeast')

figure(3)
plot(seuils_bayes, risque_bayes)
hold on
plot(seuils_eq, risque_equivalent)
grid()
xlabel("Seuil")
ylabel("Risque")
title("Risque en fonction du seuil pour Bayes et équivalent")
legend("Détecteur de Bayes", "Détecteur équivalent", 'Location', 'southeast')