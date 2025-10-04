"Ce script permet d'isoler et d'analyser les fréquences"
"dominantes dans un signal audio en utilisant la transformée de Fourier."

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

def plot_fft(signal, sample_rate):
    # Calculer la FFT
    fft_result = np.fft.fft(signal)
    # Calculer les fréquences correspondantes
    frequencies = np.fft.fftfreq(len(fft_result), 1/sample_rate)
    # Ignorer la moitié des résultats (symétrie complexe)
    positive_frequencies = frequencies[:len(frequencies)//2]
    magnitude_spectrum = np.abs(fft_result[:len(fft_result)//2])

    # Afficher le résultat
    plt.figure(figsize=(12, 6))
    plt.plot(positive_frequencies, magnitude_spectrum)
    plt.title('Spectre en fréquence')
    plt.xlabel('Fréquence (Hz)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

def analyze_audio(file_path):
    # Charger le fichier audio
    sample_rate, signal = wavfile.read(file_path)

    # Normaliser le signal audio
    signal = signal / np.max(np.abs(signal))

    # Afficher le signal temporel
    plt.figure(figsize=(12, 6))
    plt.plot(np.arange(len(signal)) / sample_rate, signal)
    plt.title('Signal audio temporel')
    plt.xlabel('Temps (s)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

    # Afficher la FFT
    plot_fft(signal, sample_rate)

# Remplacez 'chemin_du_fichier_audio.wav' par le chemin de votre fichier audio
file_path = r'.\data\croisement.wav'
analyze_audio(file_path)
