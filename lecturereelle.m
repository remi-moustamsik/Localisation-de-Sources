% Spécifiez le chemin du fichier audio
%fichier_audio = '.\data\fluteircam.wav';

% Lire le fichier audio
%[signal_audio, frequence_echantillonnage] = audioread(fichier_audio);

%Séparer les périodes où le signal peut être considéreé quasi-stationnaire
function sousListes = diviseurListe(liste, n)
%sousListes={}
% Vérifier si la longueur de la liste est un multiple de n
if mod(length(liste), n) ~= 0
    q = fix(length(liste)/n)
    r = mod(length(liste),n)
    L = liste(1:q*n)
else
    L=liste
end

% Initialiser la cellule pour stocker les sous-listes
sousListes = cell(1, fix(length(L)/n));

% Diviser la liste en sous-listes de taille n
for i = 1:length(L)/n
    debutIndex = (i - 1) * n + 1;
    finIndex = i * n;
    sousListes{i} = L(debutIndex:finIndex);
end
if mod(length(L), n) ~= 0
    sousListes{end+1} = list(end-r+1, end)
end
disp(sousListes)
end