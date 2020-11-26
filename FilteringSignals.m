%=======================
% Description
% This function comprises of low pass filter, passband filter and least
% square for filter
%
% Mfile: FilteringSignals
%
%
% Creator: Jacquline Damas
%-------------
function  [b_leastSq,Y_LowPass,Y_BandPassfiltered]=...
    FilteringSignals(OrderLow,OrderBandPass,Wcl,Wp,Ws,Rp,Rs,UnfilteredSignal,...
    LowPassFilterSelect,BandPassFiltSelect,PassbandFrequency1,...
    PassbandFrequency2,Fs,NLeastSquareOrder,Fpass,Fstop,Wpass,Wstop)
[b,a] = butter(OrderLow,Wcl,'low');               % IIR filter design
%h = fvtool(b,a);
if LowPassFilterSelect==1
    Y_LowPass = filtfilt(b,a,UnfilteredSignal);      % zero-phase filtering
else
    Y_LowPass= filter(b,a,UnfilteredSignal);         % conventional filtering
end
if BandPassFiltSelect==1
    [OrderBand,Ws] = cheb2ord(Wp,Ws,Rp,Rs);  % Gives mimimum order of filter
    [bc,ac] = cheby2(OrderBand,Rs,Ws);        % Chebyshev Type II filter
    Y_BandPassfiltered= filter(bc,ac,Y_LowPass);
else
    D = designfilt('bandpassiir', 'FilterOrder', OrderBandPass, ...
        'PassbandFrequency1',PassbandFrequency1, 'PassbandFrequency2',...
        PassbandFrequency2,...
        'SampleRate', Fs);
    Y_BandPassfiltered=filter(D,Y_LowPass);
end
% Calculate the coefficients using the FIRLS function.
b_leastSq  = firls(NLeastSquareOrder, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop]);
end

