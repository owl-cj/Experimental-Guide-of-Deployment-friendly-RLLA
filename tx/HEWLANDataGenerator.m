function [chanBW,sr,msdu,txWaveform] = HEWLANDatagenerator(numampdu,idleTime,MCS,ChannelBW,txImage,length_txImage)
%% 功能：802.11ax QoS Data MAC帧封装与PHY层波形生成
%input:
%numampdu:帧聚合数
%idleTime：每个包之间间隔时间（＞2us）
%MCS：调制编码策略
%ChannelBW：带宽
%txImage:发送载荷
%length_txImage：载荷长度
%output:
%chanBW：带宽
%sr：采样率
%msdu：MAC载荷（补零后）
%txWaveform：PHY层波形
%%
cfgSU = wlanHESUConfig;           %cfgSU: HE SU配置格式声明
cfgSU.ExtendedRange = false;      % Do not use extended-range format
cfgSU.ChannelBandwidth = ChannelBW; % Channel bandwidth
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
msduLength = 2304; % MSDU length in bytes
% numMSDUs = ceil(length(txImage)/msduLength);                  %numMSDUs:MSDU帧数计算，ceil：向上取整
numMSDUs = ceil(length_txImage/msduLength) + 1;              %后面补一帧数据长度数据帧
if numMSDUs <= numampdu                                       %numpk：WLAN包个数计算（数据会一次性发完，注意视频文件时长)
    numpk = 1;
else 
    numpk = ceil(numMSDUs/numampdu);
end
padZeros = msduLength-mod(length(txImage),msduLength);        %补零
txData = [txImage; zeros(padZeros,1)];
txDataBits = double(reshape(de2bi(txData, 8)', [], 1));       %转8bit二进制

% Divide input data stream into fragments
bitsPerOctet = 8;
data = zeros(0, 1);
pkframeind = 0;
ind2 = 0;
for ind=0:numMSDUs-1     
     pkframeind = ceil((ind+1)/ numampdu);
    % Extract image data (in octets) for each MPDU
    if ((ind2)/numampdu) == 1
        ind2 = 1;
    else
        ind2 = ind2 + 1;
    end
        if (ind == numMSDUs-1)
            msdu{pkframeind,ind2} = de2bi(length_txImage);
            else
                 msdu{pkframeind,ind2} = txData(ind*msduLength+1:msduLength*(ind+1),:);      %msdu：将txData每规定的msdu长度（2304）作为一个帧的数据部分
            end
        
end  
indnummsdu = 1;
txWaveform = zeros(0, 1);
for ind =1:numpk
    % Create MAC frame configuration object and configure sequence number
    %  cfgMAC = wlanMACFrameConfig('FrameType', 'Data', 'SequenceNumber', ind);  %帧序列号.当该MPDUAggregation属性为1 （true）时，此属性表示第一个MAC协议数据单元（MPDU）的序列号。后续MPDU的序列号以1为增量增加。
     cfgAMPDU = wlanMACFrameConfig('FrameType','QoS Data','FrameFormat','HE-SU',... 
                'MPDUAggregation',true,'MSDUAggregation',false,'SequenceNumber', indnummsdu);
     indnummsdu = indnummsdu +numampdu;
     frameBody = msdu(ind,:);
    % Generate MPDU
    [macFrames, ampduLength(ind)]= wlanMACFrame(frameBody, cfgAMPDU,cfgSU,'OutputFormat','bits');            %利用WLAN工具箱将MSDU添加相应功能单元构成MPDU(8421码 ）   
    % Concatenate PSDUs for waveform generation
    data = [data; macFrames]; %#ok<AGROW>                                %data：PHY层发送bit
    txPSDUPerUser{ind} = macFrames;
    cfgSU.APEPLength = ampduLength(ind);          % Payload length in bytes(66100maxMCS=9,A信道，过长信道的时间选择性会恶化)
    muPSDULengths(ind) = getPSDULength(cfgSU);
    % Generate HE-MU waveform with idle period
    % Initialize the scrambler with a random integer for each packet
    scramblerInitialization = randi([1 127],numMSDUs,1);
    txWaveformforpk{ind} = wlanWaveformGenerator(txPSDUPerUser{ind},cfgSU, ...
          'NumPackets',numTxPkt,'IdleTime',idleTime,'ScramblerInitialization',scramblerInitialization);
      txWaveform = [txWaveform; txWaveformforpk{1,ind}];
end
end
