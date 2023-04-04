function  MCSnew = MCSdef(snr,ChannelBWnewDec)
%% 
%功能：根据snr估计计算下一个WLAN包采用的MCS值
%输入：snr：估计的信噪比
%输出：MCSnew：计算的MCS值
%%
% SNRthreshold = [21 23 25 26 28 29 30 32 34]; %%
if ChannelBWnewDec == 4
    SNRthreshold = [10 13 15 20 23 24.5 27.5 29.5 32 34 37];
end
if ChannelBWnewDec == 2
    SNRthreshold = [11 14 16 21 23 25 28 31 33 35 37];
end
if snr>SNRthreshold(end)
    MCSnew = 11;
else if snr<SNRthreshold(1)
        MCSnew = 0;
    else
        MCSnew = find(SNRthreshold<snr, 1, 'last' );
    end
end
end
        
