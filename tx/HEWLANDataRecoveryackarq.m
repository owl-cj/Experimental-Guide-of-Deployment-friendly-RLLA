function [packetSeq,rxBitMatrix,rxnumampdu,SNRest,RXTime,pktind,searchOffset,ConstellationDiagram] = HEWLANDataRecoveryackarq(chanBW,sr,rx,idleTime)
%% 功能：WLAN HESU解包
%input:
%chanBW：带宽
%sr：采样率
%rx：接收基带波形
%output:
%packerSeq：包序列号（MAC帧序列号）
%rxBitMatrix：接收MSDU乱序
%rxnumampdu：接收计算的MAC帧聚合数
%SNRest：接收机估计的信噪比、
%%
packetSeq = [];
rxnumampdu = [];
SNRest = [];
RXTime = [];
pktind = [];
searchOffset = [];
ConstellationDiagram = [];
%接收处理，先配置接收相关参数
pilotTracking = 'Joint';
% Create an HE recovery configuration object and set the channel bandwidth
cfgRx = wlanHERecoveryConfig;
cfgRx.ChannelBandwidth = chanBW;
% The recovery configuration object is used to get the start and end
% indices of the pre-HE-SIG-B field.
ind = wlanFieldIndices(cfgRx);   %此时HE-SIG-A之前所有格式数据分布已确定，是HE packet都有的部分
% Setup plots for the example
% [spectrumAnalyzer,timeScope,ConstellationDiagram,EVMPerSubcarrier,EVMPerSymbol] = heSigRecSetupPlots(sr);
% ConstellationDiagram.Position = [400 40 1080 720];
% Minimum packet length is 10 OFDM symbols
lstfLength = double(ind.LSTF(2));
minPktLen = lstfLength*5; % Number of samples in L-STF
pktind = 1;
rxWaveLen = size(rx,1); 
summpduList = zeros(0,1);
rxBitMatrix = zeros(0,1);
packetSeq = [];
sumRXTime = 0;
%前端处理，即包检测，频率粗矫正，定时同步，频率细矫正，searchOffset超出rx索引范围时包检测失败
searchOffset = 0; % Offset from start of waveform in samples
while (searchOffset + minPktLen) <= rxWaveLen
    % Packet detection包探测，找到数据包的第一个符号，一般要提前一些，不是真正开始的符号，这么
    %做是因为方便之后频率校准的自相关操作和定时同步的精确校准
    pktOffset = wlanPacketDetect(rx,chanBW,searchOffset);  %延迟自相关算法

    % Adjust packet offset
    pktOffset = searchOffset + pktOffset;
    if isempty(pktOffset) || (pktOffset + ind.LSIG(2) > rxWaveLen)
        if pktind==1   
            disp('** No packet detected **');
            rxBitMatrix = [];
            break;
        end
    end
    %挑出传统前缀
    rxPre = rx(pktOffset+(ind.LSTF(1):ind.HESIGA(2)), :);
    %   frequency offset estimation and correction using L-STF粗频率校准
    rxLSTF = rxPre((ind.LSTF(1):ind.LSTF(2)), :);
    coarseFreqOffset = wlanCoarseCFOEstimate(rxLSTF,chanBW);
    rxPre = helperFrequencyOffset(rxPre,sr,-coarseFreqOffset);

    % Symbol timing synchronization符号定时同步（互相关算法），
    searchBufferLLTF = rxPre((ind.LSTF(1):ind.LSIG(2)),:);
    timeingoffset = wlanSymbolTimingEstimate(searchBufferLLTF,chanBW);
    pktOffset = pktOffset+timeingoffset;
    if pktOffset < 0
          disp('** No packet detected(UN) **');
          break;
    end
     rxPre = rx(pktOffset+(ind.LSTF(1):ind.HESIGA(2)), :);
     rxPre = helperFrequencyOffset(rxPre,sr,-coarseFreqOffset);

    % Fine frequency offset estimation and correction using
    % L-LTF细频率校准（SDR例程是只对前导码粗频率校准，在解得L-SIG的样本采样数后再对一个WLAN范围细频率校准，此魔改例子是对整个rx频率校准，会有处理时延等问题，有待优化）
    rxLLTF = rxPre((ind.LLTF(1):ind.LLTF(2)),:);
    fineFreqOffset = wlanFineCFOEstimate(rxLLTF,chanBW);
    rxPre = helperFrequencyOffset(rxPre,sr,-fineFreqOffset);

    % Timing synchronization complete: packet detected
%     fprintf('Packet detected at index %d\n',pktOffset + 1);

    % Display estimated carrier frequency offset
    cfoCorrection = coarseFreqOffset + fineFreqOffset; % Total CFO
%     fprintf('Estimated CFO: %5.1f Hz\n\n',cfoCorrection);



%Scale the waveform based on L-STF power (AGC)
%根据L-STF的功率归一化接收信号的功率，达到AGC效果
gain = 1./(sqrt(mean(rxLSTF.*conj(rxLSTF)))); 
rx = rx.*gain;

%包格式检测
%先对L-LTF(not include tone rotation for each 20 MHz segment)序列解调做信道估计，噪声估计，再通过这些信息从L-LTF至HE-SIG-A序列中判断包格式
rxLLTF = rxPre((ind.LLTF(1):ind.LLTF(2)),:);
lltfDemod = wlanLLTFDemodulate(rxLLTF,chanBW);
lltfChanEst = wlanLLTFChannelEstimate(lltfDemod,chanBW,11);   %信道估计(LS),使用了频率平滑
noiseVar = helperNoiseEstimate(lltfDemod);  %关于空时流的问题

% disp('Detect packet format...');
rxSIGA = rxPre((ind.LSIG(1):ind.HESIGA(2)),:);
pktFormat = wlanFormatDetect(rxSIGA,lltfChanEst,noiseVar,chanBW);
% fprintf('  %s packet detected\n\n',pktFormat);
if ~(numel(pktFormat) == 5)
      disp('** wrong target packet detected **');
      rxBitMatrix = [];
      break;
else
    regpktFormat = (pktFormat == 'HE-SU');
    if ~(regpktFormat(4)  == 1)
      disp('** wrong target packet detected **');
      rxBitMatrix = [];
      break;
    
    end
end
% Set the packet format in the recovery object and update the field indices
cfgRx.PacketFormat = pktFormat;
ind = wlanFieldIndices(cfgRx);   %由于未知是否进行HE-SIG-B压缩操作，即有可能是全带宽MU-MIMO，所以此时ind未更新
                                 %HE-SIG-B部分，若包是HE SU Ext格式，会更新HE-SIG-A长度
% 
%使用L-LTF（include tone rotation for each 20 MHz segment）做信道估计(LS)，用于在HE-LTF之前的符号做信道均衡 
lltfDemod = wlanHEDemodulate(rxLLTF,'L-LTF',chanBW);
lltfChanEst = wlanLLTFChannelEstimate(lltfDemod,chanBW,11);

%L-SIG 和 RL-SIG 的解码
% disp('Decoding L-SIG... ');
% Extract L-SIG and RL-SIG fields
rxLSIG = rxPre((ind.LSIG(1):ind.RLSIG(2)),:);

% OFDM demodulate
helsigDemod = wlanHEDemodulate(rxLSIG,'L-SIG',chanBW);

% Estimate CPE and phase correct symbols
helsigDemod = preHECommonPhaseErrorTracking(helsigDemod,lltfChanEst,'L-SIG',chanBW);

% Estimate channel on extra 4 subcarriers per subchannel and create full
% channel estimate
preheInfo = wlanHEOFDMInfo('L-SIG',chanBW);
preHEChanEst = preHEChannelEstimate(helsigDemod,lltfChanEst,preheInfo.NumSubchannels);

% Average L-SIG and RL-SIG before equalization
helsigDemod = mean(helsigDemod,2);

% Equalize data carrying subcarriers, merging 20 MHz subchannelsm（SISO默认ZF均衡，其它可选MMSE）
[eqLSIGSym,csi] = preHESymbolEqualize(helsigDemod(preheInfo.DataIndices,:,:), ...
    preHEChanEst(preheInfo.DataIndices,:,:),noiseVar,preheInfo.NumSubchannels);

% Decode L-SIG field
[~,failCheck,lsigInfo] = wlanLSIGBitRecover(eqLSIGSym,noiseVar,csi);

if failCheck
    disp(' ** L-SIG check fail **');
    rxBitMatrix = [];
end
% else
%     disp(' L-SIG check pass');
% end
% Get the length information from the recovered L-SIG bits and update the
% L-SIG length property of the recovery configuration object
lsigLength = lsigInfo.Length; %从L-SIG到RL-SIG的解码中得到了PSDU字节数
cfgRx.LSIGLength = lsigLength;

% % Measure EVM of L-SIG symbols
% EVM = comm.EVM;
% EVM.ReferenceSignalSource = 'Estimated from reference constellation';
% EVM.Normalization = 'Average constellation power';
% EVM.ReferenceConstellation = wlanReferenceSymbols('BPSK');
% rmsEVM = EVM(eqLSIGSym);
% fprintf(' L-SIG EVM: %2.1fdB\n\n',20*log10(rmsEVM/100));

%Calculate the receive time and corresponding number of samples in the
%packet
RXTime(pktind) = ceil((lsigLength + 3)/3) * 4 + 20; % In microseconds
numRxSamples = round(RXTime(pktind) * 1e-6 * sr);   % Number of samples in time
fprintf(' RXTIME of packet %d: %dus\n',pktind,RXTime(pktind));
% fprintf(' Number of samples in the packet: %d\n\n',numRxSamples);
if (numRxSamples+pktOffset)>length(rx)
    disp('**Not enough samples to decode packet **');
    break;
end

    % Apply CFO correction to the entire packet
    rx(pktOffset+(1:numRxSamples+idleTime*sr*0.75),:) = helperFrequencyOffset(...
        rx(pktOffset+(1:numRxSamples+idleTime*sr*0.75),:),sr,-cfoCorrection);
    
%在解得持续时间和采样点数后可以画出接收信号波形和频谱
% sampleOffset = max((-lstfLength + pktOffset),1); % First index to plot
% sampleSpan = numRxSamples + 2*lstfLength; % Number samples to plot
% % Plot as much of the packet (and extra samples) as we can
% plotIdx = sampleOffset:min(sampleOffset + sampleSpan,rxWaveLen);

% % Configure timeScope to display the packet
% timeScope.TimeSpan = sampleSpan/sr;
% timeScope.TimeDisplayOffset = sampleOffset/sr;
% timeScope.YLimits = [0 max(abs(rx(:)))];
% timeScope(abs(rx(plotIdx,:)));
% % release(timeScope);

% Display the spectrum of the detected packet
% spectrumAnalyzer(rx(pktOffset + (1:numRxSamples),:));
% release(spectrumAnalyzer);

% HE-SIG-A域解码
% disp('Decoding HE-SIG-A...')
rxSIGA = rx(pktOffset+(ind.HESIGA(1):ind.HESIGA(2)),:);
sigaDemod = wlanHEDemodulate(rxSIGA,'HE-SIG-A',chanBW);
hesigaDemod = preHECommonPhaseErrorTracking(sigaDemod,preHEChanEst,'HE-SIG-A',chanBW);

% Equalize data carrying subcarriers, merging 20 MHz subchannels
preheInfo = wlanHEOFDMInfo('HE-SIG-A',chanBW);
[eqSIGASym,csi] = preHESymbolEqualize(hesigaDemod(preheInfo.DataIndices,:,:), ...
                                      preHEChanEst(preheInfo.DataIndices,:,:), ...
                                      noiseVar,preheInfo.NumSubchannels);
% Recover HE-SIG-A bits
[sigaBits,failCRC] = wlanHESIGABitRecover(eqSIGASym,noiseVar,csi);

% % Perform the CRC on HE-SIG-A bits
if failCRC
    disp(' ** HE-SIG-A CRC fail **');
    rxBitMatrix = [];
    break;
end
% else
%     disp(' HE-SIG-A CRC pass');
% end

% Measure EVM of HE-SIG-A symbols
% release(EVM);
% if strcmp(pktFormat,'HE-EXT-SU')
%     % The second symbol of an HE-SIG-A field for an HE-EXT-SU packet is
%     % QBPSK.
%     EVM.ReferenceConstellation = wlanReferenceSymbols('BPSK',[0 pi/2 0 0]);
%     % Account for scaling of L-LTF for an HE-EXT-SU packet
%     rmsEVM = EVM(eqSIGASym*sqrt(2));
% else
%     EVM.ReferenceConstellation = wlanReferenceSymbols('BPSK');
%     rmsEVM = EVM(eqSIGASym);
% end
% fprintf(' HE-SIG-A EVM: %2.2fdB\n\n',20*log10(mean(rmsEVM)/100));

cfgRx = interpretHESIGABits(cfgRx,sigaBits);          
ind = wlanFieldIndices(cfgRx); % Update field indices解码HE-SIG-A后除HE-PE部分，其它部分数据分布全部确定
                                                    %可知HE-SIG-A含有较多数据包格式信息
% disp(cfgRx)

%HE-SIG-B解码
if strcmp(pktFormat,'HE-MU')
    fprintf('Decoding HE-SIG-B...\n');
    if ~cfgRx.SIGBCompression
        fprintf(' Decoding HE-SIG-B common field...\n');
        s = getSIGBLength(cfgRx);
        % Get common field symbols. The start of HE-SIG-B field is known
        rxSym = rx(pktOffset+(ind.HESIGA(2)+(1:s.NumSIGBCommonFieldSamples)),:);
        % Decode HE-SIG-B common field
        [status,cfgRx] = heSIGBCommonFieldDecode(rxSym,preHEChanEst,noiseVar,cfgRx);

        % CRC on HE-SIG-B content channels
        if strcmp(status,'Success')
            fprintf('  HE-SIG-B (common field) CRC pass\n');
        elseif strcmp(status,'ContentChannel1CRCFail')
            fprintf('  ** HE-SIG-B CRC fail for content channel-1\n **');
        elseif strcmp(status,'ContentChannel2CRCFail')
            fprintf('  ** HE-SIG-B CRC fail for content channel-2\n **');
        elseif any(strcmp(status,{'UnknownNumUsersContentChannel1','UnknownNumUsersContentChannel2'}))
            error('  ** Unknown packet length, discard packet\n **');
        else
            % Discard the packet if all HE-SIG-B content channels fail
            error('  ** HE-SIG-B CRC fail **');
        end
        % Update field indices as the number of HE-SIG-B symbols are
        % updated
        ind = wlanFieldIndices(cfgRx);
    end

    % Get complete HE-SIG-B field samples
    rxSIGB = rx(pktOffset+(ind.HESIGB(1):ind.HESIGB(2)),:);
    fprintf(' Decoding HE-SIG-B user field... \n');
    % Decode HE-SIG-B user field
    [failCRC,cfgUsers] = heSIGBUserFieldDecode(rxSIGB,preHEChanEst,noiseVar,cfgRx);

    % CRC on HE-SIG-B users
    if ~all(failCRC)
        fprintf('  HE-SIG-B (user field) CRC pass\n\n');
        numUsers = numel(cfgUsers);
    elseif all(failCRC)
        % Discard the packet if all users fail the CRC
        error('  ** HE-SIG-B CRC fail for all users **');
    else
        fprintf('  ** HE-SIG-B CRC fail for at least one user\n **');
        % Only process users with valid CRC
        numUsers = numel(cfgUsers);
    end

else % HE-SU, HE-EXT-SU
    cfgUsers = {cfgRx};
    numUsers = 1;
end

% HE-DATA 解码
cfgDataRec = trackingRecoveryConfig;
cfgDataRec.PilotTracking = pilotTracking;

% fprintf('Decoding HE-Data...\n');
for iu = 1:numUsers
    % Get recovery configuration object for each user
    user = cfgUsers{iu};  %单独的user的配置（原例程有4users)
%     if strcmp(pktFormat,'HE-MU')
%         fprintf(' Decoding User:%d, STAID:%d, RUSize:%d\n',iu,user.STAID,user.RUSize);
%     else
%         fprintf(' Decoding RUSize:%d\n',user.RUSize);
%     end

    heInfo = wlanHEOFDMInfo('HE-Data',chanBW,user.GuardInterval,[user.RUSize user.RUIndex]);

    % HE-LTF demodulation and channel estimation  
    rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
    heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',chanBW,user.GuardInterval, ...
        user.HELTFType,[user.RUSize user.RUIndex]);
    [chanEst,pilotEst] = heLTFChannelEstimate(heltfDemod,user,11);

    % Number of expected data OFDM symbols
    symLen = heInfo.FFTLength+heInfo.CPLength;
    numOFDMSym = (ind.HEData(2)-ind.HEData(1)+1)/symLen;

    % HE-Data demodulation with pilot phase and timing tracking
    % Account for extra samples when extracting data field from the packet
    % for sample rate offset tracking. Extra samples may be required if the
    % receiver clock is significantly faster than the transmitter.
    maxSRO = 120; % Parts per million
    Ne = ceil(numRxSamples*maxSRO*1e-6); % Number of extra samples
    Ne = min(Ne,rxWaveLen-numRxSamples); % Limited to length of waveform
    numRxSamplesProcess = numRxSamples+Ne;
    rxData = rx(pktOffset+(ind.HEData(1):numRxSamplesProcess),:);
    if user.RUSize==26
        % Force CPE only tracking for 26-tone RU as algorithm susceptible
        % to noise
        cfgDataRec.PilotTracking = 'CPE';
    else
        cfgDataRec.PilotTracking = pilotTracking;
    end
        [demodSym,cpe,peg] = heTrackingOFDMDemodulate(rxData,chanEst,numOFDMSym,user,cfgDataRec);

    % Estimate noise power in HE fields
    demodPilotSym = demodSym(heInfo.PilotIndices,:,:);
    nVarEst = heNoiseEstimate(demodPilotSym,pilotEst,user);

    % Equalize(MMSE)
    [eqSym,csi] = heEqualizeCombine(demodSym,chanEst,nVarEst,user);

    % Discard pilot subcarriers
    eqSymUser = eqSym(heInfo.DataIndices,:,:);
    csiData = csi(heInfo.DataIndices,:);

    % Demap and decode bits 'bp' (default) | 'layered-bp' | 'norm-min-sum' | 'offset-min-sum'
    rxPSDU = wlanHEDataBitRecover(eqSymUser,nVarEst,csiData,user,'LDPCDecodingMethod','norm-min-sum');
    Databitrate(pktind) = length(rxPSDU)/(RXTime(pktind)+idleTime*1e6);
    fprintf(' Databitrate of packet %d: %2.2f(Mbps)\n',pktind,Databitrate(pktind));
    % Deaggregate the A-MPDU
    [mpduList,AMPDUfailCRC,status] = wlanAMPDUDeaggregate(rxPSDU,wlanHESUConfig);%这里都转换成SU
     rxnumampdu = numel(mpduList);
    if strcmp(status,'Success')

    else
        fprintf('  A-MPDU deaggregation unsuccessful \n');
        rxBitMatrix = [];
        break;
    end
     summpduList = [summpduList, mpduList];
    % Decode the list of MPDUs and check the FCS for each MPDU
    for i = 1:numel(mpduList)
        [cfgheMAC,hedatapayload{i},status] = wlanMPDUDecode(mpduList{i},wlanHESUConfig,'DataFormat','octets');
        if strcmp(status,'Success')

        else
            fprintf('  FCS fail for MPDU:%d\n',i);
            rxBitMatrix = [];
            break;
        end
        
    end
    packetSeq((pktind-1)*rxnumampdu+1:pktind*rxnumampdu) = cfgheMAC.SequenceNumber-rxnumampdu+1:cfgheMAC.SequenceNumber;
    rxBitMatrix = [rxBitMatrix, hedatapayload];
    % Plot equalized constellation of the recovered HE data symbols for all
    % spatial streams per user
%     
%     hePlotEQConstellation(eqSymUser,user,ConstellationDiagram,iu,numUsers);

    % Measure EVM of HE-Data symbols
%     release(EVM);
%     EVM.ReferenceConstellation = wlanReferenceSymbols(user);
%     rmsEVM = EVM(eqSymUser(:));
%     fprintf('  HE-Data EVM:%2.2fdB\n\n',20*log10(rmsEVM/100));
%     
%     % Plot EVM per symbol of the recovered HE data symbols
%     hePlotEVMPerSymbol(eqSymUser,user,EVMPerSymbol,iu,numUsers);
%     
%     % Plot EVM per subcarrier of the recovered HE data symbols
%     hePlotEVMPerSubcarrier(eqSymUser,user,EVMPerSubcarrier,iu,numUsers);
end
% SNRset = awgnChannel.SNR;
SNRest(pktind) = 10*log10((mean(rxHELTF.*conj(rxHELTF)))/nVarEst);
% fprintf('SNRset= %2.2fdB\n',SNRset);
fprintf('SNRest=%2.2fdB\n\n ',SNRest(pktind));

searchOffset = pktOffset + ind.HEData(2);
if (length(unique(packetSeq))<length(packetSeq))  %适用于SDR
    break;
end
pktind =pktind +1;
end
if ~(isempty(rxBitMatrix))
  DataBitrate = mean(Databitrate);
  fprintf(' DataBitrate (Mbps): %2.2f\n',DataBitrate);
end
end
