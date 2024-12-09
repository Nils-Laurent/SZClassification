addpath('tftb-0.1/mfiles');
clear all; close all;
%% Crossing chirps:

rng(0);
N = 2^10;
Nchirp = N;
tmin = round((N-Nchirp)/2);
tmax = tmin + Nchirp;
x = zeros(N,1);
x(tmin+1:tmax) = (fmlin(Nchirp,0.15,0.35) + fmlin(Nchirp,0.35,0.15)).*tukeywin(Nchirp,0.25);

% Parameters for the STFT.
Nfft = 2*N;
fmax = 0.5; % Max. norm. frequency to compute the STFT.
[w,T] = roundgauss(Nfft,1e-6); % Round Gaussian window.

% Noise realization:
SNRin = 20;
noise = randn(size(x));
xnoise = sigmerge(x,noise,SNRin);
[F,~,~] = tfrstft(xnoise,1:N,Nfft,w,0);
F = F(1:floor(Nfft*fmax),:);
F = flipud(F);
S = abs(F).^2;
% imagesc(-abs(F).^0.3); hold on;

% Find original zeros and triangulation
ceros = find_zeros_stft(S);
TRI = delaunay(ceros);

% Keep zeros within margins:
margin_row = 10; margin_col = 10;
invalid_ceros = zeros(length(ceros),1);
invalid_ceros(ceros(:,1)<margin_row | ceros(:,1)>(size(S,1)-margin_row))=1;
invalid_ceros(ceros(:,2)<margin_col | ceros(:,2)>(size(S,2)-margin_col))=1;
invalid_ceros = logical(invalid_ceros);
valid_ceros = ~invalid_ceros;
% number_of_valid_ceros = sum(valid_ceros);

% Noise assisted approach:
disp('Computing histogram...');
lims = 0;
aux_plane = zeros(size(S));
M = 16; % Number of noise realizations.
SNRalg = SNRin;
% SNRest = estimate_SNR(xnoise,S,N,Nfft,w);

parfor j  = 1:M
    noise_alg = randn(N,1);
    xnoise_alg = sigmerge(xnoise,noise_alg,SNRalg);
    [S_alg,~,~] = tfrsp(xnoise_alg,1:N,Nfft,w,0);
    S_alg = S_alg(1:floor(Nfft*fmax),:);
    S_alg = flipud(S_alg);
    [~, Qz] = find_zeros_stft(S_alg);
    aux_plane = aux_plane + Qz;
end

hist2d = aux_plane(lims+1:end-lims,lims+1:end-lims);
selected_hist = zeros(size(hist2d));

% Eliminate some spourious acumulations on the border of the plane.
selected_hist(5+1:end-5,5+1:end-5) = hist2d(5+1:end-5,5+1:end-5);
selected_hist = selected_hist/sum(selected_hist(:));

disp('Finished.');

%% Show the histogram.
% figure()
% imagesc(selected_hist)

%% Compute a concentration measure on patches around zeros and a threshold.
dT = 2; ceil(T/4);
% Compute a KD Tree for easily searching balls of neighbors:
[row,col] = ind2sub(size(selected_hist),1:numel(selected_hist));
idx = [row.' col.'];
Mdl = KDTreeSearcher(idx);

% figure()
% Compute the entropy on patches around zeros.
for i = 1:size(ceros,1)
    % Squared Patches:
    %     coordx = max([ceros(i,2)-dT,1]):min([ceros(i,2)+dT,size(selected_hist,2)]);
    %     coordy = max([ceros(i,1)-dT,1]):min([ceros(i,1)+dT,size(selected_hist,1)]);
    %     patch = selected_hist(coordy,coordx);
    
    % Circular Patches:
    patch_ind = rangesearch(Mdl,ceros(i,:),dT);
    patch = selected_hist(patch_ind{1});
    
    %     Uncomment this to see the circular patches:
    %     aux = zeros(size(selected_hist));
    %     aux(patch_ind{1}) = 1;
    %     imagesc(aux); hold on;
    %     plot(ceros(:,2),ceros(:,1),'o','Color','r','MarkerFaceColor','r','MarkerSize',4); hold on;
    %     clf;
    
    Hceros(i) = max(patch(:));
%     Hceros(i) = sum(patch(:));
    
end

%%
% Compute threshold and keep zeros above it
Hthr = quantile(Hceros,0.98);
% Hthr = mean(Hceros);
ind_selected_zeros = Hceros>Hthr;
selected_zeros = ceros(ind_selected_zeros,:);

%
figure()
% subplot(1,2,1)
imagesc(-abs(F).^0.3); hold on;
colormap bone;
plot(ceros(:,2),ceros(:,1),'o','Color','w','MarkerFaceColor','w','MarkerSize',4);
plot(selected_zeros(:,2), selected_zeros(:,1),'o','Color','r','MarkerFaceColor','r','MarkerSize',4);
xticks([])
yticks([])
xlabel('time')
ylabel('frequency')
title(sprintf('Spectrogram and zeros (%2.0fdB)',SNRalg))


figure()
% subplot(1,2,2)
imagesc(hist2d.^0.3);
xticks([])
yticks([])
xlabel('time')
ylabel('frequency')
title('Selected Zeros')
