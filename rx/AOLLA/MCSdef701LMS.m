function  [MCSnew,PERList,deltaThresholdList] = MCSdef701LMS(snr,ChannelBWnewDec,ctlinfoList,PERList,deltaThresholdList)
%% 
%功能：根据snr估计计算下一个WLAN包采用的MCS值
%输入：snr：估计的信噪比
%输出：MCSnew：计算的MCS值
%%
global nt;
if isempty(nt)
    nt=0;
end
global u;
global deltaThreshold;
if nt<= 9999
    PERest = 10*log10(0.092);
    disp('10%');
else
    PERest = 10*log10(0.05);
     disp('5%');
end
if length(ctlinfoList) <= 8
    PER = length(find(ctlinfoList==2))/length(ctlinfoList);
end
if length(ctlinfoList) > 8
    PER = length(find(ctlinfoList(end-7:end)==2))/8;
end

disp(['PER'  num2str(PER) ] );
PERList = [PERList PER];
if PER ==0
    PER = -20;
%     PER = -16;
else
    PER = 10*log10(PER);
end
    e = PER - PERest;
    if e >0
        e = 1.75*e;
    end
    if e<0
        e = 0.65*e;
    end
if isempty(u)
    u = 0.1;
end
if isempty(deltaThreshold)
    deltaThreshold = 0;
end
deltaThreshold = deltaThreshold + u * e;
u = 0.002 * (u + 3 * abs(e));
disp(['deltaThreshold ' num2str(deltaThreshold)]);
deltaThresholdList = [deltaThresholdList deltaThreshold];
% SNRthreshold = [21 23 25 26 28 29 30 32 34]; %%
if ChannelBWnewDec == 4
    SNRthreshold = [10 13 15 20 23 24.5 27.5 29.5 32 34 37] + deltaThreshold;
%     SNRthreshold = [12 16 17 21 24 26 29.5 31 35 37 40] + deltaThreshold;
end
if ChannelBWnewDec == 2
    SNRthreshold = [11 14 16 21 23 25 28 31 33 35 37] + deltaThreshold;
%     SNRthreshold = [12 16 17 21 24 26 29.5 31 35 37 40] + deltaThreshold;
end
if snr>SNRthreshold(end)
    MCSnew = 11;
else if snr<SNRthreshold(1)
        MCSnew = 0;
    else
        MCSnew = find(SNRthreshold<snr, 1, 'last' );
    end
end
nt = nt+1;
disp(['nt=' num2str(nt)]);
end
        
