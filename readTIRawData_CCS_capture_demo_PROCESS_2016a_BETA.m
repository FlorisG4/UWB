% Parse raw TI radar data captured via CCS
load('CorrectionCoefficients_colours.mat');
RadarID = input('Select Radar ID: 1 (green), 2 (yellow), 3 (blue), 4 (red): ');
switch RadarID
    case 1
        CorCoef = CorCoef_Green;
    case 2
        CorCoef = CorCoef_Yellow;
    case 3
        CorCoef = CorCoef_Blue;
    case 4
        CorCoef = CorCoef_Red;
    otherwise
        disp('Invalid Radar ID');
        return
end

NTxProc = input('Select Number of Tx cahnnels for AOA processing (1 or 2): ');
if NTxProc>2
    display('Invalid number of Tx antennas')
    return
end
%% Load radar settings from powershell script
[fileName,path] = uigetfile('*.ps1','Select the script file containing the radar config'); % select data file to load
fid = fopen([path fileName]);
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% channel config
%channelCfgStr = string(c{1,1}(~cellfun(@isempty,strfind(c{1,1},'channelCfg'))));
channelCfgStr = char(c{1,1}(~cellfun(@isempty,strfind(c{1,1},'channelCfg'))));
channelCfg = str2double((regexp(channelCfgStr,'(\d+,)*\d+(\.\d*)?', 'match')));

rxChannelEn = channelCfg(1);
txChannelEn = channelCfg(2);
nRxEn = length(strfind(dec2bin(rxChannelEn),'1'));
nTxEn = length(strfind(dec2bin(txChannelEn),'1'));

% profile config
profileCfgStr = char(c{1,1}(~cellfun(@isempty,strfind(c{1,1},'profileCfg'))));
profileCfg = str2double((regexp(profileCfgStr,'(\d+,)*\d+(\.\d*)?', 'match')));

startFr = profileCfg(2);            % GHz
idleTime = profileCfg(3);           % microseconds 
adcStartTime = profileCfg(4);       % microseconds
rampEndTime = profileCfg(5);        % microseconds
freqSlopeConst = profileCfg(8);     % MHz/microseconds
txStartTime = profileCfg(9);        % microseconds
numAdcSamples = profileCfg(10);
digOutSampleRate =  profileCfg(11); % kHz

% chirp config
chirpCfgStr = char(c{1,1}(~cellfun(@isempty,strfind(c{1,1},'chirpCfg'))));
nChirps = size(chirpCfgStr,1);
for n = 1:nChirps
    chirpCfg(n,:) = str2double((regexp(chirpCfgStr(n,:),'(\d+,)*\d+(\.\d*)?', 'match')));
end

% frame config
frameCfgStr = char(c{1,1}(~cellfun(@isempty,strfind(c{1,1},'frameCfg'))));
frameCfg = str2double((regexp(frameCfgStr,'(\d+,)*\d+(\.\d*)?', 'match')));
nTxOn = abs(frameCfg(2)-frameCfg(1))+1;
nLoops = frameCfg(3);
nFrames = frameCfg(4);
framePeriodicity = frameCfg(5); % ms

%% Load TI RAW data file
NSamp = numAdcSamples; % samples per chirp
NChirps = nLoops; % chirps per loop per Tx channel
NTx = nTxOn;
NRx = nRxEn;
NChan = NTx*NRx;
NFrames = nFrames;

Fs = digOutSampleRate*1e3; % Hz
dt = 1/Fs;
t = 0:dt:(NSamp-1)*dt; % seconds

[fileName,path] = uigetfile('*.raw','Select the data file'); % select data file to load
fid = fopen([path fileName],'r');
adcRaw = fread(fid,'int16');
fclose(fid);

adcRaw = adcRaw(1:2:end)+1i*adcRaw(2:2:end); % interleaved IQ

% reshape data array
adcRaw2 = reshape(adcRaw,[NRx,NTx*NSamp,NChirps,NFrames]);
if NTx==2
    adcRaw2 = [adcRaw2(:,1:NSamp,:,:);adcRaw2(:,NSamp+1:end,:,:)];
end

% figure;
% subplot(4,1,1);plot(t*1e6,squeeze(real(radar_data(1,:,1))),t*1e6,squeeze(imag(radar_data(1,:,1))))
% legend('IRx1','QRx1');xlabel('t (\mus)');title('Beat signals 1.PostProc')
% subplot(4,1,2);plot(t*1e6,squeeze(real(radar_data(2,:,1))),t*1e6,squeeze(imag(radar_data(2,:,1))))
% legend('IRx2','QRx2');xlabel('t (\mus)');
% subplot(4,1,3);plot(t*1e6,squeeze(real(radar_data(3,:,1))),t*1e6,squeeze(imag(radar_data(3,:,1))))
% legend('IRx3','QRx3');xlabel('t (\mus)');
% subplot(4,1,4);plot(t*1e6,squeeze(real(radar_data(4,:,1))),t*1e6,squeeze(imag(radar_data(4,:,1))))
% legend('IRx4','QRx4');xlabel('t (\mus)');

%% Range profiles (one frame)
NFFTR = 1024;
winR = repmat(ones(NChan,1)*hann(NSamp).',[1,1,NChirps]);
RP = fft(winR.*squeeze(adcRaw2(:,:,:,1)),NFFTR,2);

fb = linspace(0,Fs,NFFTR);
S = freqSlopeConst*1e12; % Slope [Hz/s]
Range = 3e8*fb/(2*S);
% figure;
% subplot(4,1,1);plot(Range,squeeze(20*log10(abs(RP(1,:,1)))));
% legend('Rx1');xlabel('range (m)');title('Range profiles 1.PostProc')
% subplot(4,1,2);plot(Range,squeeze(20*log10(abs(RP(2,:,1)))));
% legend('Rx2');xlabel('range (m)');
% subplot(4,1,3);plot(Range,squeeze(20*log10(abs(RP(3,:,1)))));
% legend('Rx3');xlabel('range (m)');
% subplot(4,1,4);plot(Range,squeeze(20*log10(abs(RP(4,:,1)))));
% legend('Rx4');xlabel('range (m)');

% %% Range Doppler (one channel)
% NFFTD = 256;
% winRD = hann(NSamp)*hann(NChirps)';
% RD = fft2(winRD.*squeeze(adcRaw2(1,:,:,1)),NFFTR,NFFTD);
% 
% figure;imagesc([],Range,fftshift(20*log10(abs(RD)/max(abs(RD(:)))),2));axis xy;
% colormap jet;caxis([-40 0])
%% Range / AOA
NFFTA = 128;
AOA = asind( linspace(-1, 1, NFFTA) );
NChan = NTxProc*4;
NoCorCoef = ones(1,NChan);
winAOA = ones(NChan,1)*ones(1,NFFTR);

% Plotting
RAOA = fft(winAOA.*(squeeze(RP(1:NChan,:,1)).*(CorCoef(1:NChan).'*ones(1,NFFTR))),NFFTA);
figure; hp = pcolor(AOA,Range,fftshift(20*log10(abs(RAOA)/max(abs(RAOA(:))))',2));
set(hp,'EdgeColor','none'); caxis([-40 0])
colormap jet;xlabel('angle (deg)');ylabel('range (m)');
title('Range/Angle map')
ylim([1 3])

% % Angle cut at selected range
uiwait(msgbox('Click on figure to select 1D cut','modal'));
[x,y]=ginput(1);

% MUSIC start
objects = 3;
distance = y;                                                              % Distance of object, estimated from graph
A = winAOA.*(squeeze(RP(1:NChan,:,1)).*(CorCoef(1:NChan).'*ones(1,NFFTR))); % Calibrating the range data
r_bin = round(distance/max(Range)*NFFTR);                                      
[U,S,V] = svd(A(:, r_bin-4:r_bin+4));                                       % Doing the SVD
elementPos = [-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75];
direction = zeros(length(AOA), 1);                                          % Initialize vector to store direction data


for i=1:length(AOA)
    ang = [AOA(i);0];
    sv = steervec(elementPos,ang);
    direction(i) = 1/sum(abs(U(:, 1+objects:end)'*sv).^2);
end
plot(AOA, db(direction))
grid()
xlabel('Direction (°)')
ylabel('Magnitude')
title(sprintf('MUSIC on distance %0.2f m', distance));
% Music END

% Regular Angle Cut start

[minv,mini]=min(abs(Range-y));
figure;plot(AOA,fftshift(20*log10(abs(RAOA(:,mini))/max(abs(RAOA(:))))',2))
grid;
xlabel('angle (deg)');ylabel('mag (dB)');title(['Angle cut at ' num2str(y,3) ' m'])
% 
% [~,maxi] = max(abs(RAOA(:)));
% [~,J] = ind2sub(size(RAOA),maxi);
% figure;plot(AOA,fftshift(20*log10(abs(RAOA(:,J))/max(abs(RAOA(:))))',2))
% grid;
% xlabel('angle (deg)');ylabel('mag (dB)');title(['Angle cut at ' num2str(Range(J),3) ' m (maximum)'])

% Regular Angle Cut end
