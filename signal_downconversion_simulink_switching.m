clear all; close all; clc;
%% --- Description---------------------------------------
%
%
%Mfile: Signal_downconversion.m
%
%
% Creator: Jacqueline Damas

N=1;
for k=1:N
    %% Initialize system parameters
    %Fc = 2e9;                                           % Carrier frequency in GHz
    Fc = 2.0000119997e9;
    %Fc = 2.00001197e9;
    %Fc = 2.0000125e9;
    %Fc = 999999046.33;
    %Fc=1.999998093e9;
    Fs =20e9;                                           % Signal sampling frequency in GHz
    Ts = 1/Fs;                                          % Sampling time interval in s
    Rs= 1e9;                                            % Symbol rate
    sps = Fs/Rs;                                        % Samples per symbol
    T = 1/Rs;                                           % Symbol time interval in s
    rolloff = 0.5;
    M = 4;                                              % Modulation order (QPSK)
    span = 10;
    % L = Fs/Rs;                                        % Oversampling factor
    InitVal=10.47e5;                                    % remove unwanted samples
    %% Read csvread files
    if k==1
        data = csvread('1Gbps_Data_20Gbps.csv');
        %data = data(InitVal:end,:);
    elseif k==2
        data = csvread('1Gbps_Data_20Gbps_to_MZM.csv');
        data = data(InitVal:end,:);
    elseif k==3
        data = csvread('1Gbps_Data_40Gbps_detected.csv');
        data = data(InitVal:end,:);
        Fs =40e9;
    elseif k==4
        data = csvread('1Gbps_Data_40Gbps_to_MZM.csv');
        data = data(InitVal:end,:);
        Fs =40e9;
    elseif k==5
        data = csvread('3Gbps_Data_20Gbps_detected.csv');
        data = data(InitVal:end,:);
    elseif k==6
        data = csvread('3Gbps_Data_40Gbps_detected.csv');
        data = data(InitVal:end,:);
        Fs =40e9;
    elseif k==7
        data = csvread('3Gbps_Data_40Gbps_MZM.csv');
        data = data(InitVal:end,:);
        Fs =40e9;
    else
        data = csvread('3Gbps_Data_40Gbps_to_MZM.csv');
        data = data(InitVal:end,:);
        Fs =40e9;
    end
    t = data(:,1);
    amp = data(:,2);
    
    %%
    UnfilteredSignal=data(:,2);
    %--Low pass filter
    OrderLow=15;
    Wcl=3.1e9/(0.5*Fs);
    
    %% Band pass parameters
    %-----Standard bandpass param----------
    OrderBandPass=20;
    PassbandFrequency1=250000;
    PassbandFrequency2=3.0e9;
    
    if (k==1||k==3) % noise signals-- detected
        %----Cheby2 Parameters----------
        Wcl=3.1e9/(0.5*Fs);
        Wch1=1e9;
        Offset1=0.7e9;
        Offset2=0.7e9;
        Wch2=3.1e9;
        Factr_Fs=0.5;
        Wp = [Wch1 Wch2]/(Factr_Fs*Fs); %pass band between 0.7e5 to 3.2e9
        Ws = [(Wch1-Offset1) (Wch2+Offset2)]/(Factr_Fs*Fs); %0.6e9Hz wide on both sides of the passband
        Rp = 3;   %3 dB of ripple
        Rs = 18;  %  27dB attenuation and above is worse
    else
        Wch1=30000000;
        Wch2=3.2e9;
        Factr_Fs=0.5;
        Wp = [Wch1 Wch2]/(Factr_Fs*Fs); %pass band between 0.7e5 to 3.2e9
        Ws = [25000000 3.8e9]/(Factr_Fs*Fs); %0.6e9Hz wide on both sides of the passband
        Rp = 3;   %3 dB of ripple
        Rs = 7;  % 40 dB attenuation
    end
    %% Design parameters LPF with least square
    NLeastSquareOrder     = 25;  % Order
    Fpass = 1.1e9;               % Passband Frequency
    Fstop = 1.112e9;             % Stopband Frequency
    Wpass = 1;                   % Passband Weight
    Wstop = 1;                   % Stopband Weight
    
    %% FilterTypeSelect
    
    BandPassFiltSelect=1; % 1= Cheby2 otherwise other bandpass
    LowPassFilterSelect=2; % 1 Zero_Phase fitering otherwise conventional
    
    %% filtering
    
    [b_leastSq,Y_LowPass,Y_BandPassfiltered]=...
        FilteringSignals(OrderLow,OrderBandPass,Wcl,Wp,Ws,Rp,Rs,UnfilteredSignal,...
        LowPassFilterSelect,BandPassFiltSelect,PassbandFrequency1,...
        PassbandFrequency2,Fs,NLeastSquareOrder,Fpass,Fstop,Wpass,Wstop);
   
    %% Demodulation
    %Phi_1=27.5;
    %Phi_2=27.5;
    %Phi_1=25.3;
    %Phi_2=25.3;
    Phi_1=118.9;
    Phi_2=118.9;
    IConv = sin(2*pi*Fc*t+Phi_1);
    QConv = -cos(2*pi*Fc*t+Phi_2);
    r_I = IConv .* Y_BandPassfiltered;
    r_Q = QConv .* Y_BandPassfiltered;
    
            %% spectrum of the filtered signal
    spectrum = Y_BandPassfiltered;
    %filter = filtfilt(b_leastSq, 1, Y_BandPassfiltered);
    %data = filter;
    P=length(spectrum);
    t = 0:1:P-1;
    phi = zeros(1,P);
    R=1:P;
    N = [0 : P-1];
    N_2 = ceil(N/2);
    Q = -P/2:P/2-1;
    Mag = abs(fft(spectrum));
    Mag_d = 20*log10(Mag);
    Mag_dB = 20*log10(Mag/max(Mag)); % Normalised magnitude by the maximum
    %X = Mag/P;
    %S = Fs * length(R-1); %number of samples
    Freq_Hz = N*Fs/P;
    P_2 = ceil(P/2);
    
    
    
    
%     frameSize = P;
%     filterOversample = 6;
%     frequencyOffsetHz = 1e8;
%     % Shift signal in frequency
%     t = 0:1/Fs:(frameSize*filterOversample-1)/Fs;
%     freqShift = exp(1i.*2*pi*frequencyOffsetHz*t.');
%     offsetData = Y_BandPassfiltered.*freqShift;
    
    
    
    %% filtering after demodulation
    w_I = filtfilt(b_leastSq, 1, r_I);
    w_Q = filtfilt(b_leastSq, 1, r_Q);
    
    %% eye diagram with filtering
    Factr=4.5;
    B_I_eye=reshape(w_I(1:floor(length(w_I)/(sps*Factr))*(sps*Factr)),[(sps*Factr),floor(length(w_I)/(sps*Factr))]);
    B_Q_eye=reshape(w_Q(1:floor(length(w_Q)/(sps*Factr))*(sps*Factr)),[(sps*Factr),floor(length(w_Q)/(sps*Factr))]);
    nrsample=sps*10;
    u_I = downsample(w_I, nrsample, nrsample/2);
    u_Q = downsample(w_Q, nrsample, nrsample/2);
    
    
    
    %sim_I = B_I_eye;
    %sim_Q = B_Q_eye;
    % %sim_T = sim_TTT(:,1:100);
    %sim_T = sim_TTT(1:1000,1:2);
    
    
    sim_I = w_I;
    sim_Q = w_Q;
    siz_I = size(sim_I,1);
    siz_Q = size(sim_Q,1);
    tt = 5e-11;
    n_I = siz_I*tt;
    n_Q = siz_Q*tt;
    t_I = linspace(0,n_I,siz_I);
    t_Q = linspace(0,n_Q,siz_Q);
    sim_II = timeseries(w_I(:,1),t_I);
    sim_QQ = timeseries(w_Q(:,1),t_Q);
    
    %sim_II = timeseries(sim_I(:,1),t_I);
    %sim_QQ = timeseries(sim_Q(:,1),t_Q);
    
    %sim_II = timeseries(sim_I(1:90,1:2));
    %sim_QQ = timeseries(sim_Q(1:90,1:2));
    %sim_II =timeseries(sim_I(1:90,2),sim_I(1:90,1));
    %sim_QQ =timeseries(sim_Q(1:90,2),sim_Q(1:90,1));
    
    %sim_T = sim_TTT[(1:100 1:100)(1:100 1:100)];
    %sim_T = ([1:100] sim_TTT(:,1:100));
    
    
%     figure;
%     subplot(411);plot([1:(sps*Factr)],B_I_eye(:,1:150));
%     if k==1
%         title('1Gbps-Data-20Gbps-detected')
%     elseif k==2
%         title('1Gbps-Data-20Gbps-to-MZM')
%     elseif k==3
%         title('1Gbps-Data-40Gbps-detected')
%     else
%         title('1Gbps-Data-40Gbps-to-MZM-')
%     end
%     grid on
%     hold on
%     subplot(412);plot([1:(sps*Factr)],B_Q_eye(:,1:150));
%     grid on
%     legend('Eye')
%     % %constellation
%     %subplot(313);plot(u_I,u_Q,'.');
%     subplot(413);plot(B_I_eye(1,1:1500),B_Q_eye(1,1:1500),'.');
%     grid on
%     legend('constellation')
%     subplot(414);plot(Freq_Hz(1:P_2), Mag_dB(1:P_2));
%     xlim([0 4*10^9])
%     grid on
    
end

%plot (w_I(1:1000,1)) %%plotting I downconverted with the respect to number of samples
%plot (w_Q(1:1000,1))  %%plotting Q downconverted with the respect to number of samples

%plot(abs(fft(data)));
%plot(abs(fftshift(data)));


% xdft = fft(data);
% xdft(1:length(data)/2+1);
% [Y,I] = max(abs(xdft));
% freq = 0:Fs/length(data):Fs/2;
% fprintf('Maximum occurs at %3.2f Hz\n.',freq(I))


%Y=fft(y);
%X=fftshift(Y);
%f=linspace(-Fs/2,+Fs/2,length(X));
%plot(f,20*log10(abs(X)));
%</code>

