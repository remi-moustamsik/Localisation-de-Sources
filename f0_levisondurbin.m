fichier_audio = ".\data\croisement.wav";
fe=32000;f0=440;f1=500;
xx=cos(2*pi*f0/fe*[1:1280]+2*pi*rand(1,1))+cos(2*pi*f1/fe*[1:1280]+2*pi*rand(1,1))
[signal_audio, frequence_echantillonnage] = audioread(fichier_audio);

N = 2^16;

n = 0.004*frequence_echantillonnage;

nt = floor(length(signal_audio)/n);

sub_lists = subdivide_signal(signal_audio, nt);
%test pour les signaux enregistrés
% n[aar, sig2, refl, fdsp, dsp] = mylevisondurbin(sub_lists{200}', 500, fe);
[liste1 , liste2] = plot_frqz(signal_audio, floor(nt/20), fe)
%test pour 2 cosinus de fréquences proches
%[aar, sig2, refl, fdsp, dsp] = mylevisondurbin(xx , 500, fe);

function [list_freq1 , list_freq2] = plot_frqz(signal_audio, k, fe)

    sub_lists = subdivide_signal(signal_audio, k)
    
    list_freq1 = cell(k , 1)
    list_freq2 = cell(k , 1)
    list_ampl1 = cell(k , 1)
    list_ampl2 = cell(k , 1)

    for i = 1:length(sub_lists)
        [aar, sig2, refl, fdsp, dsp] = mylevisondurbin(sub_lists{i}', 500, fe)
        [f01 , f02 , a01 , a02] = find_f0_2(sub_lists{i}')
        list_ampl1{i} = a01
        list_ampl2{i} = a02
        % Utiliser la fonction findpeaks pour trouver les maxima locaux
        %[pics, indices] = findpeaks(dsp);
        % Trier les maxima locaux par ordre décroissant
        %double_liste = [pics', indices'];

        % Trier la double liste en fonction de la première colonne (pics)
        %sorted_double_liste = sortrows(double_liste, 1);

        % Séparer les vecteurs triés
        %sorted_pics = sorted_double_liste(:, 1);
        %sorted_indices = sorted_double_liste(:, 2);

        %créer les listes des fréquences
        if i == 1
            list_freq1{i} = f01
            list_freq2{i} = f02
        if i ~= 1 
            if abs(f01-f02)<0.1f01
                if abs(list_ampl1{i-1}-a01)<=abs(list_ampl2{i-1}-a01)
                    list_ampl1{i} = a01
                    list_ampl2{i} = a02
                else
                    list_ampl1{i} = a02
                    list_ampl2{i} = a01
                end
            end
        else
            if abs(list_freq1{i-1}-f01)<abs(list_freq2{i-1}-f01)
                list_freq1{i} = f01
                list_freq2{i} = f02
            else
                list_freq2{i} = f01
                list_freq1{i} = f02
            end
        end
    end
    % Afficher le graphique avec les deux courbes
    plot(cell2mat(list_freq1)); % Courbe 1 en bleu avec une épaisseur de ligne de 2
    hold on; % Permet de superposer les graphiques
    plot(cell2mat(list_freq2)); % Courbe 2 en rouge en pointillés avec une épaisseur de ligne de 2
    end
% Ajouter des labels, une légende et un titre
xlabel('temps');
ylabel('fréquence');
title('suivie des fréquences par rapport au temps');
end

function [f01 , f02 , a1 , a2] = find_f0_2(trame)
    [~, ~, ~, freqdsp, dsp] = mylevinsondurbin(trame', 100, 3200);
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

function [aa, sigma2, ref, ff, mydsp] = mylevisondurbin (xx, pp, fe)

acf = xcorr(xx, pp+1, 'biased'); %% autocorr\'elation
acf(1:pp+1) = [];                %% on enl\`eve la partie n\'egative
acf(1) = real(acf(1));           %% Levinson-Durbin requiert c(1)==conj(c(1))

ref = zeros(pp,1);
gg = -acf(2)/acf(1);
aa = [ gg ];
sigma2 = real( ( 1 - gg*conj(gg)) * acf(1) ); %% real : enl\`eve une \'eventuelle
                                              %%   partie imaginaire r\'esiduelle
ref(1) = gg;
for tt = 2 : pp
  gg = -(acf(tt+1) + aa * acf(tt:-1:2)') / sigma2;
  aa = [ aa + gg*conj(aa(tt-1:-1:1)), gg ];
  sigma2 = sigma2 * ( 1 - real(gg*conj(gg)) );
  ref(tt) = gg;
end;
aa = [1, aa];

%%% densit\'e spectrale de puissance
interm2=-j*2*pi/fe*[1:pp];
df=0.9765625;      %%% la dsp est calcul\'ee tous les df Hz
ff=-fe/2:df:fe/2;

interm3=interm2'*ff;
interm=1.+aa(2:pp+1)*exp(interm3);
mydsp = sigma2./(interm.*conj(interm));

%figure(1);
%clf;
%grid on;
%hold on;
%semilogy(ff,mydsp,'linewidth',2);
%xlabel('frequency (in Hz)','fontsize',20);
%ylabel('magnitude','fontsize',20);
%hold off;
%drawnow;

end

function subdivided_lists = subdivide_signal(signal_audio, k)
    
    signal_length = length(signal_audio);
    
    subdivision_size = floor(signal_length / k);
    
    subdivided_lists = cell(k, 1);
    
    for i = 1:k
        start_index = (i - 1) * subdivision_size + 1;
        end_index = min(i * subdivision_size, signal_length);
        
        subdivided_lists{i} = signal_audio(start_index:end_index);
    end
end

function f0 = find_f0(i,fft_list,fe, m)
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

function smoothed_list = sliding_window_average(input_list, window_size)
    
    smoothed_list = NaN(size(input_list));
    
    for i = 1:length(input_list)
        start_index = max(1, i - floor(window_size/2));
        end_index = min(length(input_list), i + floor(window_size/2));
        
        window_values = input_list(start_index:end_index);
        smoothed_list(i) = mean(window_values, 'omitnan');
    end
end