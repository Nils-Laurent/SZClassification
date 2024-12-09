clear all
addpath('tftb-0.1/mfiles');
close all
clc

%% Signals
N = 2^8;
tmin = round(sqrt(N));
tmax = N-tmin;
Nchirp = round(tmax-tmin);
tchirp = (0:Nchirp-1);

semilla=rng();
tone1 = cos(2*pi*(0.27*tchirp));
tone2 = cos(2*pi*(0.23*tchirp));
xmid = (tone1 + tone2).'.*tukeywin(Nchirp,0.25);
x = zeros(N,1);
x(tmin+1:tmax) = xmid;

% Parameters for the STFT.
fmax = 0.5; % Max. norm. frequency to compute the STFT.
Nfft = 2*N;
[w,T] = roundgauss(Nfft,1e-6); % Round Gaussian window.

% [S,~,~] = tfrsp(x,1:N,Nfft,w);
% figure()
% imagesc(flipud(S(1:Nfft/2,:)));
% axis square

SNRin = 30:-0.1:-10;
noise = randn(size(x));
signal = sigmerge(x,noise,SNRin(1));

%% Original zeros and triangulation
[mask,sr,F,TRI,original_ceros,TRIselected,w,T,TRI_EDGES] = delaunay_domains(signal);
% Enumerate the zeros
zero_id = 1:length(original_ceros);
S = abs(F).^2;

% figure()
% imagesc(-abs(F).^0.3); hold on;
% colormap bone
% % plot(original_ceros(:,2),original_ceros(:,1),'+',...
% %     'Color','r','MarkerFaceColor','r','MarkerSize',5); hold on;
% triplot(TRI,original_ceros(:,2),original_ceros(:,1),'Color','b')
% for i = 1:length(original_ceros)
%     text(original_ceros(i,2),original_ceros(i,1),sprintf('%d',i));
% end

%%
figure()
imagesc(abs(F)); hold on;
plot(original_ceros(:,2),original_ceros(:,1),'o',...
    'Color','r','MarkerFaceColor','r','MarkerSize',5); hold on;

% TRI = delaunay(original_ceros(:,1),original_ceros(:,2));
% [vx_original,vy_original] = voronoi(original_ceros(:,1),original_ceros(:,2),TRI);
% plot(vy_original,vx_original,'b--','LineWidth',0.2); hold off;
axis square
FRAME(1) = getframe;

for i = 2:length(SNRin)
    % Compute the new zeros for a fix realization and different variance 
    xnoise = sigmerge(x,noise,SNRin(i));
    [F,~,~] = tfrstft(xnoise,1:N,Nfft,w,0);
    F = F(1:floor(Nfft*fmax),:);
    F = flipud(F);
    [ceros, Qz] = find_zeros_stft(abs(F));    
    imagesc(abs(F)); hold on;
    
    % Plot original ceros, the voronoi cells and the new zeros
    plot(original_ceros(:,2),original_ceros(:,1),'o',...
        'Color','r','MarkerFaceColor','r','MarkerSize',5);
%     plot(vy_original,vx_original,'b--','LineWidth',0.2)
    plot(ceros(:,2),ceros(:,1),'o','MarkerSize',5,'MarkerFaceColor','auto');
    axis square

%     snr_frame = 20*log10(norm(x)/norm(new_noise(:,i)+noise));
    snr_frame = SNRin(i);
    text(5,5,sprintf('frame #%d - snr = %2.1f dB ',i,snr_frame),'FontSize',15);
    text(original_ceros(:,2)+0.1,original_ceros(:,1),string(1:length(original_ceros)),'FontSize',10);
    FRAME(i) = getframe;
    hold off;
end

save('frame.mat','FRAME');
figure()
movie(FRAME)


