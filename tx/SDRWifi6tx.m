%%
clear all;
% Check that WLAN Toolbox is installed, and that there is a valid
% license
if isempty(ver('wlan')) % Check for WLAN Toolbox install
    error('Please install WLAN Toolbox to run this example.');
elseif ~license('test', 'WLAN_System_Toolbox') % Check that a valid license is present
    error( ...
        'A valid license for WLAN Toolbox is required to run this example.');
end
% % Setup Spectrum viewer
% spectrumScope = dsp.SpectrumAnalyzer( ...
%     'SpectrumType', 'Power density', ...
%     'SpectralAverages', 10, ...
%     'YLimits', [-130 -40], ...
%     'Title', 'Received Baseband WLAN Signal Spectrum', ...
%     'YLabel', 'Power spectral density');
%%
%  Initialize SDR device
deviceNameSDR = 'Pluto'; % Set SDR Device
radio = sdrdev(deviceNameSDR);           % Create SDR device object
% %% Input an image file and convert to binary stream
% fileTx = 'E:\毕业设计\图片素材\test2\images\001.jpg';            % Image file name
% [imsize,txImage,length_txImage] = imgframeinput(fileTx);
%% 得到视频压缩后的bin数据 并初始化分块
% workingdir = 'C:\Users\chenjie\Desktop\matlabcode\yournamejpg';
% imageNamesTx = dir(fullfile(workingdir,'*.xlsx'));
% imageNamesTx = {imageNamesTx.name}';
% lengthfile = xlsread('C:\Users\chenjie\Desktop\matlabcode\yournamejpg\lengthfile.xlsx');
% for i = 1:numel(imageNamesTx)-1
%     fileTx = fullfile(workingdir,imageNamesTx{i});
%     txImagefile(:,i) = xlsread(fileTx);
% end
length_txcode = 7610368*100;
txcode = ceil(256*rand(length_txcode,1)) -  1;
% fid = fopen('BBB1080.bin','rb');
% [txcode,length_txcode] = fread(fid,'uint8');
% fclose(fid);

% 将发送视频序列以20个WLAN包为单位分段
numAPk = 1;
numampdu = 1;                    %一个WLAN包最多聚合的msdu数，与信道的时间选择性有关
numampdunew = 12; 
timesTx = ceil(length_txcode/(2304*numAPk*numampdu));
for i = 1:timesTx-1  
    txCodeAPk{i} =  txcode(2304*numAPk*numampdu*(i-1)+1:2304*numAPk*numampdu*i);
end
txCodeAPk{timesTx} = txcode(2304*numAPk*numampdu*(timesTx-1)+1:end);

%% Input an image file and convert to binary stream (frameind：发送图片帧指示）
frameind = 1;
% imsize = [720 1280 3];
MCS = 4;
ChannelBW = 'CBW40';
ChannelBWnew = 'CBW40';

while (frameind<=timesTx)
%     fileTx = fullfile(workingdir,imageNamesTx{frameind});
%     txImage = xlsread(fileTx);
%      txImage = txImagefile(:,frameind);
%     length_txImage = lengthfile(frameind);
txImage = [];
numampdu = numampdunew;                    %一个WLAN包最多聚合的msdu数，与信道的时间选择性有关
if numampdu == 1
     txImage = txCodeAPk{frameind} ;
else
     for ampduind = 1:numampdu
        txImage = [txImage; txCodeAPk{frameind+ampduind-1}] ;
     end
end
length_txImage = length(txImage);
%% 将DATA按照WLAN标准形成发射波形

ChannelBW = ChannelBWnew;
idleTime = 2e-6; % Idle time before and after each packet
[chanBW,sr,msdu,txWaveform] = HEWLANDatagenerator(numampdu,idleTime,MCS,ChannelBW,txImage,length_txImage);

%% The sdrTransmitter uses the |transmitRepeat| functionality to transmit
% the baseband WLAN waveform in a loop from the DDR memory on the PlutoSDR.
sdrTransmitter = sdrtx(deviceNameSDR); % Transmitter properties
sdrTransmitter.RadioID = 'usb:0';

% Resample the transmit waveform at 30 MHz
fs = sr; % Transmit sample rate in MHz
osf = 1.5;                     % Oversampling factor

sdrTransmitter.BasebandSampleRate = fs* osf; 
sdrTransmitter.CenterFrequency = 3.184e9;  % Channel 5
sdrTransmitter.ShowAdvancedProperties = true;
sdrTransmitter.Gain = 0;

% Resample transmit waveform 
txWaveform  = resample(txWaveform,fs*osf,fs);

fprintf('\nGenerating WLAN transmit waveform:\n')

% Scale the normalized signal to avoid saturation of RF stages
powerScaleFactor = 0.8;
txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor);

% Transmit RF waveform
 sdrTransmitter.transmitRepeat(txWaveform);
 disp(['CBW' num2str(ChannelBWnew) num2str(numampdunew)]);
disp(frameind);
pause(8); 
 release(sdrTransmitter); 
%% 接收ack/arq包并选择下一个发送的WLAN包的信息
rxBitMatrixackarq = [];
BasebandSampleRate = fs*osf;   %SDR采样率
 CenterFrequency = 1.784e9;     %中心频率
% *Receiver Setup*
sdrReceiverackarq = sdrrx(deviceNameSDR);
sdrReceiverackarq.RadioID = 'usb:0';
sdrReceiverackarq.BasebandSampleRate = BasebandSampleRate;
sdrReceiverackarq.CenterFrequency = CenterFrequency;
sdrReceiverackarq.GainSource = 'AGC Slow Attack';
% sdrReceiver.Gain = 16;
sdrReceiverackarq.OutputDataType = 'double';

captureLength = 1e6;%60MHz采样率采12M个点，即采样持续0.2s
while isempty(rxBitMatrixackarq)
fprintf('\nStarting a new ack/arq RF capture.\n')    
        burstCapturesackarq = capture(sdrReceiverackarq, captureLength, 'Samples');
% Downsample the received signal
rxackarq = resample(burstCapturesackarq,fs,fs*osf);
% processing WLAN packet
[~,rxBitMatrixackarq,~,~,~,~,~,~] = HEWLANDataRecoveryackarq(chanBW,sr,rxackarq,idleTime);
if ~(isempty(rxBitMatrixackarq))
if ~(isempty(rxBitMatrixackarq{1}))
%     spectrumScope.SampleRate = BasebandSampleRate;
%     spectrumScope(burstCapturesackarq);
%     release(spectrumScope);

      rxinfoackarq = hex2dec(rxBitMatrixackarq{1}{1});
      ackarqinfo = char(rxinfoackarq(1:3));
      MCSnext = rxinfoackarq(4);
      numampdunew = rxinfoackarq(5);
      if numampdunew > 12
          numampdunew = 12;
      end
      ChannelBWnewDec = rxinfoackarq(6);
      if ChannelBWnewDec == 2
           ChannelBWnew = 'CBW20';
      else 
           ChannelBWnew = 'CBW40';
      end
else
%     if ~(isempty(rxBitMatrixackarq{2}))
%       rxinfoackarq = hex2dec(rxBitMatrixackarq{2}{1});
%       ackarqinfo = char(rxinfoackarq(1:3));
%       MCSnext = rxinfoackarq(4);
%       numampdunew = rxinfoackarq(5);
%       ChannelBWnewDec = rxinfoackarq(6);
%       if ChannelBWnewDec == 2
%            ChannelBWnew = 'CBW20';
%       else 
%            ChannelBWnew = 'CBW40';
%       end
%     end
%     if ((isempty(rxBitMatrixackarq{1}))&&(isempty(rxBitMatrixackarq{2})))
        MCSnext = 3;
        numampdunew = 12;
        ChannelBWnewDec = 4;
        ackarqinfo = 'arq';        
    end
end
end

release(sdrReceiverackarq);
%% ack/arq决策
if (ackarqinfo(2) == 'c')&&(ackarqinfo(3) == 'k')
    frameind = frameind +1;
    disp('ack');
    
else
    if (ackarqinfo(2) == 'r')&&(ackarqinfo(3) == 'q')
        disp('arq');
    end
end
MCS = MCSnext;

fprintf('MCSnext = %d \n',MCS);
end