 function  [txWaveform,MCSnew,msdu,PERList,deltaThresholdList,numampdunew] = ackarqtx1104QLearning(ctlinfo,snr,numampdunew,ChannelBWnewDec,chanBW,ctlinfoList,PERList,deltaThresholdList,MCSbuffer,beta1,beta2)
feedbackind = [1,3,4,5,6];
asciiack = abs('ack');
asciiarq = abs('arq');
if (length(ctlinfoList) <= 5 || (ctlinfo == 0))
    MCSnew = MCSdef(snr,ChannelBWnewDec);
    if ctlinfo == 0
        MCSnew = MCSbuffer;
    end
end
if (length(ctlinfoList) > 5 && ctlinfo ~= 0)  
    [MCSnew,PERList,deltaThresholdList] = MCSdef1104QLearning(snr,ChannelBWnewDec,ctlinfoList,PERList,deltaThresholdList,beta1,beta2);
end
    
% if (length(ctlinfoList) <= 10 || ctlinfo == 0)  %前10个包不采用卡尔曼，10这个数要调参，这是该方案的的固有缺点：
%    MCSnew = MCSdef(snr,ChannelBWnewDec);        %初始数据数量太少，导致在前面的数据包输出产生激进的MCS选择策略，造成误包，
% end                                             %数据太多，造成卡尔曼无法起到效果。也指示了测试时选择良好的初始信道条件
% %% 观测矩阵数值选取
% if length(ctlinfoList) > 10 
%   MCSnew = MCSdef630Kalman(snr,ChannelBWnewDec,ctlinfoList,snravg); 
% end
if ctlinfo == 1
    msdu = [asciiack,MCSnew];
    disp('ack');
      fprintf('new MCS = %d\n',MCSnew);
else
    if ~(MCSnew == 0)
    MCSnew = MCSnew ;
    end
    msdu = [asciiarq,MCSnew];
      fprintf('new MCS = %d\n',MCSnew);
    disp('arq');
end
if ((MCSnew == 3)||(MCSnew == 2))&&(ChannelBWnewDec == 2)
    if numampdunew >6
        numampdunew = 6;
    end
end
if (MCSnew == 1)&&(ChannelBWnewDec == 2)
    if numampdunew >3
        numampdunew = 3;
    end
end
if (MCSnew == 0)&&(ChannelBWnewDec == 2)
    if numampdunew >1
        numampdunew = 1;
    end
end    
if (MCSnew == 1)&&(ChannelBWnewDec == 4)
    if numampdunew >6
        numampdunew = 6;
    end
end
if (MCSnew == 0)&&(ChannelBWnewDec == 4)
    if numampdunew >3
        numampdunew = 3;
    end
end
  msdu = [msdu,numampdunew,ChannelBWnewDec];
%%
MCS = 0;
idleTime = 2e-6;
cfgSU = wlanHESUConfig;           %cfgSU: HE SU配置格式声明
cfgSU.ExtendedRange = false;      % Do not use extended-range format
cfgSU.ChannelBandwidth = chanBW; % Channel bandwidth
cfgSU.MCS = MCS;                   % Modulation and coding scheme
cfgSU.ChannelCoding = 'LDPC';     % Channel coding
cfgSU.NumSpaceTimeStreams = 1;    % Number of space-time streams
cfgSU.NumTransmitAntennas = 1;    % Number of transmit antennas
numpk = 1;                        %需要发送的WLAN包个数

                                  %通过wlanWaveformGenerator得到HE SU包
% Specify waveform parameters
numTxPkt = 1; % Number of transmitted packets
chanBW = cfgSU.ChannelBandwidth;
sr = wlanSampleRate(cfgSU);

% *Fragment transmit data*
     cfgAMPDU = wlanMACFrameConfig('FrameType','QoS Data','FrameFormat','HE-SU',... 
                'MPDUAggregation',false,'MSDUAggregation',false);
     frameBody = msdu;
    % Generate MPDU
    [macFrames, ampduLength]= wlanMACFrame(frameBody, cfgAMPDU,cfgSU,'OutputFormat','bits');   %利用WLAN工具箱将MSDU添加相应功能单元构成MPDU(8421码 ）   
    % Concatenate PSDUs for waveform generation
    txPSDUPerUser = macFrames;
    cfgSU.APEPLength = ampduLength;          % Payload length in bytes(66100maxMCS=9,A信道，过长信道的时间选择性会恶化)
    muPSDULengths = getPSDULength(cfgSU);
    % Generate HE-MU waveform with idle period
    % Initialize the scrambler with a random integer for each packet
    scramblerInitialization = randi([1 127],1,1);
    txWaveform = wlanWaveformGenerator(txPSDUPerUser,cfgSU, ...
          'NumPackets',numTxPkt,'IdleTime',idleTime,'ScramblerInitialization',scramblerInitialization);

end