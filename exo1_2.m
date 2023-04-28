clc; clear; close all

load("data.mat")

%fonction = exp(-(1/2)*((X-v)/sigma).^2) .* exp((1/2)*(X/sigma).^2); %detecteur de bayes
numerator = zeros(N, n);
denominator = zeros(N, n);

for i = 1:N
    for j = 1:n
        numerator(i, j) = sum(log(1/sqrt(2*pi*sigma^2)) - 0.5*((Z(i,j)-X(i))/sigma).^2);
        denominator(i, j) = sum(log(1/sqrt(2*pi*sigma^2)) - 0.5*((Z(i,j))/sigma).^2);
    end
end

numerator = sum(log(numerator), 2);
denominator = sum(log(denominator), 2);

fonction = numerator ./ denominator;
lambda_0 = (pi0*(c10-c00))/(pi1*(c01-c11)); %seuil detecteur de bayes
seuil = 2*sigma^2*log(lambda_0) + sum(X.^2); %seuil detecteur equivalent

d1 = 1;
d0 = 0;

%Detecteur de Bayes
delta1 = zeros(N, n);

for i = 1:N
    compteur_Bayes = 0;
    for j = 1:n
        if (fonction(1, j) > lambda_0)
            delta1(i, j) = d1;
            compteur_Bayes = compteur_Bayes + 1;
        else
            delta1(i, j) = d0;
        end
    end
end

%Detecteur equivalent
delta2 = zeros(N, n);

for i = 1:N
    compteur_equivalent = 0;
    for j = 1:N
        if (Z(i,j) > seuil)
            delta2(i, j) = d1;
            compteur_equivalent = compteur_equivalent + 1;
        else
            delta2(i, j) = d0;
        end
    end
end

for i = 1:N
    nb_d1_detecte_bayes = 0;
    nb_fa_bayes = 0;
    nb_d1_detecte_equivalent = 0;
    nb_fa_equivalent = 0;
    
    for j = 1:n
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

delta_eq = (X > fonction);
nb_d1_detecte = sum(delta_eq .* vraies_detection);
nb_fa = sum(delta_eq .* ~vraies_detection);
Pd = nb_d1_detecte ./ sum(vraies_detection);
Pfa = nb_fa ./ sum(~vraies_detection);

figure(1)
plot(X(1,:))
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