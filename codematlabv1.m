% Spécifiez le chemin du fichier audio
fichier_audio = '.\data\fluteircam.wav';

% Lire le fichier audio
[signal_audio, frequence_echantillonnage] = audioread(fichier_audio);

% Calculer la transformée de Fourier du signal audio sur une portion
% quasi-constante
transformee_fourier = fft(signal_audio);

% Obtenir la fréquence associée à chaque composante de la transformée
frequence = (0:length(transformee_fourier)-1) * (frequence_echantillonnage / length(transformee_fourier));

% Tracer le spectre en magnitude
figure;
plot(frequence, abs(transformee_fourier));
title('Spectre en magnitude de la transformée de Fourier');
xlabel('Fréquence (Hz)');
ylabel('Magnitude');
grid on;
