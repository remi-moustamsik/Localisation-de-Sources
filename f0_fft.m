fichier_audio = '.\data\croisement.wav';

[signal_audio, frequence_echantillonnage] = audioread(fichier_audio);

N = 2^16;

n = 0.004*frequence_echantillonnage;

nt = floor(length(signal_audio)/n);

sub_lists = subdivide_signal(signal_audio, nt);

ffts = cell(nt, 1);

engs = (1:nt);

for i = 1:nt
    ffts{i} = abs(fft(sub_lists{i}, N));
end

for i = 1:nt
    engs(i) = autocorr_0(sub_lists{i});
end

%figure;
%plot((1:nt), engs);
%title('Energie par trame');
%xlabel('Trame (numero)');
%ylabel('energie');
%grid on;
    

%Obtenir la fréquence associée à chaque composante de la transformée
%frequence = (0:length(transformee_fourier)-1) * (frequence_echantillonnage / length(transformee_fourier));

% Tracer le spectre en magnitude
%figure;
%plot(frequence, abs(transformee_fourier));
%title('Spectre en magnitude de la transformée de Fourier');
%xlabel('Fréquence (Hz)');
%ylabel('Magnitude');
%grid on;



%test = find_f0(15, ffts,frequence_echantillonnage,N);

function subdivided_lists = subdivide_signal(signal_audio, k) %Découper l'enregistrement en trames
    
    signal_length = length(signal_audio);
    
    subdivision_size = floor(signal_length / k);
    
    subdivided_lists = cell(k, 1);
    
    for i = 1:k
        start_index = (i - 1) * subdivision_size + 1;
        end_index = min(i * subdivision_size, signal_length);
        
        subdivided_lists{i} = signal_audio(start_index:end_index);
    end
end

function f0 = find_f0(i,fft_list,fe, m) %Déterminer la fréquence fondamentale d'une trame
    transformee_fourier = fft_list{i};
    frequence = (0:(m-1)) * (fe / m);

    tf = transformee_fourier(1:(length(transformee_fourier)/2));
    freq = frequence(1:(length(frequence)/2));

    [~, index] = max(tf);
    
    f0 = freq(index);
    
end

function eng = autocorr_0(trame)
    autocorr = xcorr(trame);
    eng = max(autocorr);
end 

function smoothed_list = sliding_window_average(input_list, window_size) %lisser le signal
    
    smoothed_list = NaN(size(input_list));
    
    for i = 1:length(input_list)
        start_index = max(1, i - floor(window_size/2));
        end_index = min(length(input_list), i + floor(window_size/2));
        
        window_values = input_list(start_index:end_index);
        smoothed_list(i) = mean(window_values, 'omitnan');
    end
end

function [f01 , f02 , a1 , a2] = find_f0_2(trame) 
    [~, ~, ~, freqdsp, dsp] = mylevinsondurbin(trame', 100, frequence_echantillonnage);
    tf = sliding_window_average(dsp(1:(length(dsp)/2)), 5);
    freq = freqdsp(1:(length(freqdsp)/2));
    
    global_max = max(tf);
    
    modified_list = tf;
    modified_list(modified_list < 0.5 * global_max) = -1;
    
    i0 = 1;
    flag = 0;
    lastpos = 0;

    while flag == 0
        val = modified_list(i0);
        if val > 0
            if lastpos == 0
                lastpos = 1;
            end
        end

        if val < 0
            if lastpos == 1
                flag = 1;
            end
        end
        i0 = i0 + 1;
    end

    [~, index] = max(modified_list(1:i0));

    f01=freq(index)
    
    a01=tf(index);

    i1=i0

    flag = 0;
    lastpos = 0;

    while flag == 0
        val = modified_list(i0);
        if val > 0
            if lastpos == 0
                lastpos = 1;
            end
        end

        if val < 0
            if lastpos == 1
                flag = 1;
            end
        end
        i0 = i0 + 1;
    end

    [~, index] = max(modified_list(i1:i0));

    f02=freq(index)
    
    a02=tf(index);

end
