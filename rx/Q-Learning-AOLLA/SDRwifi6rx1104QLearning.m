%%
clear all;
%%

% Check that WLAN Toolbox is installed, and that there is a valid
% license
if isempty(ver('wlan')) % Check for WLAN Toolbox install
    error('Please install WLAN Toolbox to run this example.');
elseif ~license('test', 'WLAN_System_Toolbox') % Check that a valid license is present
    error( ...
        'A valid license for WLAN Toolbox is required to run this example.');
end
% Setup Spectrum viewer
% spectrumScope = dsp.SpectrumAnalyzer( ...
%     'SpectrumType', 'Power density', ...
%     'SpectralAverages', 10, ...
%     'YLimits', [-130 -40], ...
%     'Title', 'Received Baseband WLAN Signal Spectrum', ...
%     'YLabel', 'Power spectral density');

% Setup the constellation diagram viewer for equalized WLAN symbols
% constellation = comm.ConstellationDiagram(...
%     'Title', 'Equalized WLAN Symbols', ...
%     'ShowReferenceConstellation', false);
%% Receiver Design
PERList = [];
deltaThresholdList = [];
frameind = 1;
ackind = 1;
timesTx = 200;
idleTime = 2e-6;
rxdata = [];
%  Initialize SDR device
deviceNameSDR = 'Pluto'; % Set SDR Device
radio = sdrdev(deviceNameSDR);           % Create SDR device object
chanBW = 'CBW40' ;                %基带带宽
sr = 40000000;                    %采样率
% imsize = [720 1280 3];             %图片分辨率
 fs = sr;                                   %采样率
 osf = 1.5;                     % Oversampling factor
 BasebandSampleRate = fs*osf;   %SDR采样率
 CenterFrequency = 3.184e9;     %中心频率
% *Receiver Setup*
sdrReceiver = sdrrx(deviceNameSDR);
sdrReceiver.RadioID = 'usb:0';
sdrReceiver.BasebandSampleRate = BasebandSampleRate;
sdrReceiver.CenterFrequency = CenterFrequency;
sdrReceiver.GainSource = 'AGC Slow Attack';
% sdrReceiver.Gain = 16;
sdrReceiver.OutputDataType = 'double';

captureLength = 1.2e7;%60MHz采样率采12M个点，即采样持续0.2s
% spectrumScope.SampleRate = sdrReceiver.BasebandSampleRate;
SNRperPk = zeros(100,1);
rmsEVMsymDiff = [];
rmsEVMscDiff = [];
ChannelBWnewDec = 4;
feedbackState = zeros(4,200);  %行一ack/arq，行二MCS，行三MSDU，行四BW
feedbackState(1,1) = 1;
feedbackState(2,1) = 4;
feedbackState(3,1) = 12;
feedbackState(4,1) = 4;
ctlinfoList = [];
feedbackQ = zeros(5,1);
Nt = 0;
MCSbuffer = 4;
WQLearing = 25;
PacketDuration(:,1) = [1348;708;484;388;276;228;196;196;164;164;148;132];               %不同帧聚合等级,不同MCS下包长统计
PacketDuration(:,2) = [3908;1988;1348;1028;708;548;484;452;388;356;324;292];
PacketDuration(:,3) = [0;3908;2628;1988;1348;1028;916;836;708;644;580;516];
PacketDuration(:,4) = [0;0;3908;2948;1988;1508;1348;1220;1028;932;836;756];
PacketDuration(:,5) = [0;0;5188;3908;2628;1988;1764;1604;1348;1220;1092;980];
        k = 1;
actionSpace = [-0.2 0 0.2];
stateSpace = [0.2 0.4 0.6 0.8 1.2 1.4 1.6 1.8 2.0 2.2 2.4];
load Qtable1114.mat
QtableOld = Qtable;
QtableSNR = 25;
vQtable1 = Qtable(2:4,2:5); 
meanQtable1 = mean(vQtable1(:));
vQtable2 = Qtable(2:4,6:12); 
meanQtable2 = mean(vQtable2(:));
QvalueNew1 = meanQtable1;
QvalueNew2 = meanQtable2;
Qtable5 = QtableOld;
% Qtable = zeros(length(actionSpace)+1,length(stateSpace)+1)+10000;
% Qtable(2:4,1)= actionSpace';
% Qtable(1,2:12)= stateSpace';
aQtable1 = Qtable(2:4,2:5);
[action1lo,beta1lo] = find(Qtable == max(aQtable1(:)));
action1(1) = Qtable(action1lo,1);
beta1(1) = Qtable(1,beta1lo);
aQtable2 = Qtable(2:4,6:12);
[action2lo,beta2lo] = find(Qtable == max(aQtable2(:)));
action2(1) = Qtable(action2lo,1);
beta2(1) = Qtable(1,beta2lo);
beta1new(1) = beta1(1);
beta2new(1) = beta2(1);
% action1(1) = 0;
% action2(1) = 0;
% beta1(1) = 0.4;
% beta2(1) = 1.8;
% beta1new(1) = 0.4;
% beta2new(1) = 1.8;
gama = 0.9;
epsilon = 1;
epsilonTarget = 0.01;
kWarmUp = 1;
alpha1(1) = 0.1;
alpha2(1) = 0.1;
alpha1min = 0.1;
alpha2min = 0.1;
beta2up = 0;
beta2down = 0;
beta1up = 0; 
beta1down = 0 ;
while(Nt<=timesTx)
fprintf('\nStarting a new RF capture.\n')
if ChannelBWnewDec == 2
    sr = 20000000;
else 
    sr = 40000000;
end
fs = sr;                       %采样率
BasebandSampleRate = fs*osf;   %SDR采样率
sdrReceiver.BasebandSampleRate = BasebandSampleRate;    
burstCaptures = capture(sdrReceiver, captureLength, 'Samples');
release(sdrReceiver);
% *Receiver Processing* 
% spectrumScope(burstCaptures);
%  release(spectrumScope);
% Downsample the received signal
rx = resample(burstCaptures,fs,fs*osf);

%% WLAN解包 MCS,帧聚合，BW的决策
ctlinfo = 2;
[packetSeq,rxBitMatrix,rxnumampdu,SNRest,RXTime,pktind,searchOffset,~,ctlinfo,rmsEVMsym,rmsEVMsc] = HEWLANDataRecovery(chanBW,sr,rx,idleTime);
if  ctlinfo ~= 0
    ctlinfoList = [ctlinfoList;ctlinfo];
end
disp(['recive msdu:' num2str(rxnumampdu)]);
disp(['total tx time(ms) :' num2str(0.001*sum(RXTime))]);
SNRack = mean(SNRest);
disp(['SNRack(dB) :' num2str(SNRack)]);
numampdunew = 12;
ChannelBWnewDec = 4;
if ctlinfo ~= 0  %因为同步问题导致的丢包，数据无意义,选择丢弃它
    SNRperPk(frameind) = mean(SNRest);
    %帧聚合决策：每个MSDU段SNR差距大于2dB就减少一级聚合度
     SNRtime=[];
     SNRSwap = ceil(length(rmsEVMsym)/rxnumampdu);
     if isempty(rmsEVMsym)       % 尺寸帧未通过CRC
            SNRtime = 0;
            rmsEVMsc = zeros(468,1);
     else
        for i = 1:ceil(length(rmsEVMsym)/SNRSwap)-1
            SNRtime(i) = mean(rmsEVMsym(1+(i-1)*SNRSwap:i*SNRSwap));
        end
     end
            if isempty(i)               %聚合数为1时
                SNRtime = mean(rmsEVMsym);
            else
                if ~(isempty(rmsEVMsym))
                 SNRtime(i+1) = mean(rmsEVMsym(1+i*SNRSwap:end));
                end
            end
  
            if max(SNRtime)-min(SNRtime)>5
                if (rxnumampdu == 3 || rxnumampdu == 1) 
                    numampdunew =1;
                else
                    numampdunew = rxnumampdu-3;
                end
        
            end
          if ctlinfo == 1
            if max(SNRtime)-min(SNRtime)<=5
                if rxnumampdu == 1
                    numampdunew = 3;
                else if rxnumampdu < 12
                        numampdunew = rxnumampdu + 3;
                        if rxnumampdu >= 12
                            numampdunew = 12;
                        end
                    end
                end
            end
          end
%             else
%                 if rxnumampdu == 1
%                     numampdunew = 3;
%                 else if rxnumampdu ~= 12
%                         numampdunew = rxnumampdu + 3;
%                     end
%               
%             end
        disp(['SNR vs time dB :' num2str(max(SNRtime)-min(SNRtime))]);
        disp(['newMSDU :' num2str(numampdunew)]);
    %带宽决策：带内外子载波SNR均值差距大于10dB缩小带宽
        if chanBW(4) == '4'
          BWdiff = abs(mean([rmsEVMsc(1:117);rmsEVMsc(352:468)])-mean(rmsEVMsc(118:351)));
            if abs(mean([rmsEVMsc(1:117);rmsEVMsc(352:468)])-mean(rmsEVMsc(118:351))) >= 12
            ChannelBWnewDec = 2;
%             chanBW = 'CBW20' ;                %基带带宽
 
            
            end
       disp([ num2str(BWdiff) 'dB:']);
        end
        if chanBW(4) == '2'
            if  ctlinfo == 1
                ChannelBWnewDec = 4;
            end
        end


    disp(['newBW(10MHz) :' num2str(ChannelBWnewDec)]);   
    rmsEVMsymDiff = [rmsEVMsymDiff max(SNRtime)-min(SNRtime)];
    rmsEVMscDiff = [rmsEVMscDiff BWdiff];
end

%% MAC帧按序号排列，解码出MSDU数据
if ctlinfo == 1 || ctlinfo ==2
Nt = Nt+1;
disp(Nt);
end
if ctlinfo == 1
% rxData = imgframeoutput(packetSeq,rxBitMatrix,rxnumampdu,pktind,searchOffset,sr);
% %     % Plot received image
% %         figure; 
% %     imshow(receivedImage);
% %     title(sprintf('Received Image'));
% rxdata = [rxdata rxData'];
 frameind = frameind + 1;
 disp(frameind-1);
end
%% 发送ACK/ARQ 数据包
% [ackarqtxWaveform,MCSnew,ackarqmsdu] = ackarqtx(ctlinfo,SNRack,numampdunew,ChannelBWnewDec,chanBW);
% [ackarqtxWaveform,MCSnew,ackarqmsdu,PERList,deltaThresholdList] = ackarqtx701LMS(ctlinfo,SNRack,numampdunew,ChannelBWnewDec,chanBW,ctlinfoList,PERList,deltaThresholdList);
[ackarqtxWaveform,MCSnew,ackarqmsdu,PERList,deltaThresholdList,numampdunew] = ackarqtx1104QLearning(ctlinfo,SNRack,numampdunew,ChannelBWnewDec,chanBW,ctlinfoList,PERList,deltaThresholdList,MCSbuffer,beta1new(k),beta2new(k));
if ctlinfo ~= 0
    MCSbuffer = MCSnew;
end
ackind = ackind+1;
if ChannelBWnewDec == 4
    chanBW = 'CBW40';
else
    chanBW = 'CBW20';
end
     feedbackState(1,ackind) = ctlinfo;
     feedbackState(2,ackind) = MCSnew;
     feedbackState(3,ackind) = numampdunew;
     feedbackState(4,ackind) = ChannelBWnewDec;
     feedbackState(5,ackind-1) = SNRack;%0项是未同步丢包，计算时不予考虑
     
     if ctlinfo ~= 0
        feedbackQ(:,Nt) =  feedbackState(:,ackind-1);
     end
    if length(feedbackQ(1,:)) == 10+WQLearing*k

        for kk = 1:WQLearing
            rewardState(:,kk) = feedbackQ(:,Nt-WQLearing+kk);
        end
        reward(k) = fRewardCal (rewardState);
        if k > kWarmUp
            if k <= 2 * kWarmUp
            epsilon = ((epsilonTarget-1)/kWarmUp)*k+1-(epsilonTarget-1);
            else
                epsilon = epsilonTarget;
            end
        end
        if k > 1
           Qtable = fQtableCal (reward(k),Qtable,action1(k),action2(k),beta1(k),beta2(k),beta1new(k),beta2new(k),gama,alpha1(k),alpha2(k));
        end
        k = k + 1;
        disp(['k = '  num2str(k)]);
        beta1(k) = beta1new(k-1);
        beta2(k) = beta2new(k-1);
%        [QtableBuffer,Qtable0,Qtable1,Qtable2,Qtable3,Qtable4,Qtable6,Qtable7,Qtable8,Qtable9,Qtable10,Qtable5]= fQtableBuffer(Qtableold,Qtable,rewardState,QvalueNew1,QvalueNew2,beta1new(k-1),beta2new(k-1),lo1,lo2,lo11,lo3,Qtable0,Qtable1,Qtable2,Qtable3,Qtable4,Qtable6,Qtable7,Qtable8,Qtable9,Qtable10,Qtable5);
       [action1(k),action2(k)] = fActionUpdate(epsilon,Qtable,beta1new(k-1),beta2new(k-1));
       beta1new(k) = action1(k) + beta1new(k-1);
       beta1new(k) = round(beta1new(k)*10)/10;
       beta2new(k) = action2(k) + beta2new(k-1);
       beta2new(k) = round(beta2new(k)*10)/10;
       if beta1new(k)<0.2
           beta1new(k) = 0.2;
           beta1down = beta1down + 1;
       end
       if beta1new(k)>0.8
           beta1new(k) = 0.8;
           beta1up = beta1up + 1;
       end
       if beta2new(k)<1.2
           beta2new(k) = 1.2;
           beta2down = beta2down + 1;
       end
       if beta2new(k)>2.4
           beta2new(k) = 2.4;
           beta2up = beta2up + 1;
       end
       
       %学习率alpha递减
%        [alpha1(k),alpha2(k)] = fAlpha(epsilon,epsilonTarget,action1(k),action2(k),beta1(k),beta2(k));
%        [alpha1(k),alpha2(k)] = fAlpha(epsilon,epsilonTarget,alpha1min,alpha2min);
       alpha1(k) = 0.1;
       alpha2(k) = 0.1;
    end
% the baseband WLAN waveform in a loop from the DDR memory on the PlutoSDR.
sdrTransmitterackarq = sdrtx(deviceNameSDR); % Transmitter properties
sdrTransmitterackarq.RadioID = 'usb:0';

% Resample the transmit waveform at 30 MHz
fs = sr; % Transmit sample rate in MHz
osf = 1.5;                     % Oversampling factor

sdrTransmitterackarq.BasebandSampleRate = fs*osf; 
sdrTransmitterackarq.CenterFrequency = 1.784e9;  % Channel 5
sdrTransmitterackarq.ShowAdvancedProperties = true;
sdrTransmitterackarq.Gain = 0;

% Resample transmit waveform 
ackarqtxWaveform  = resample(ackarqtxWaveform,fs*osf,fs);

fprintf('\nGenerating WLAN transmit waveform:\n')

% Scale the normalized signal to avoid saturation of RF stages
powerScaleFactor = 0.8;
ackarqtxWaveform = ackarqtxWaveform.*(1/max(abs(ackarqtxWaveform))*powerScaleFactor);

% Transmit RF waveform
sdrTransmitterackarq.transmitRepeat(ackarqtxWaveform);
if ctlinfo == 2
    wait_time = 7;
else
    wait_time =5;
end
pause(wait_time);
release(sdrTransmitterackarq);

end
%% 保存接收到的视频编码文件
% % fid2 = fopen('D:\videocodec\HM\HM\workspace\hevcOutput.bin','wb');
% fid2 = fopen('C:\Users\24390\Desktop\hevcOutput.bin','wb');
% count2 = fwrite(fid2,rxdata');
