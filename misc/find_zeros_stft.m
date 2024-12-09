function [ceros, Qz] = find_zeros_stft(S)
% Find zeros as local minima, i.e. in a grid of 3x3, of the spectrogram S
% Input:
% - S: Spectrogram (squared magnitude of the STFT)
% Output:
% -ceros: A Nx2 array, where N is the number of zeros.
% - Qz: A matrix with the same size as S, with a "1" where a cero is 
% located and "0" otherwise.

% Juan Manuel Miramont-Taurel
%--------------------------------------------------------------------------
N1=size(S,1);
M1=size(S,2);
k=1;
saltProx=0;
Qz = zeros(size(S));
% Inf=max(max(S));

S=[Inf*ones(N1,1) S Inf*ones(N1,1)];
S=[Inf*ones(1,M1+2); S; Inf*ones(1,M1+2)];

N=size(S,1);
M=size(S,2);

for i = 2:N-1
    for j=2:M-1
        
        if saltProx==1
            saltProx=0;
        else
            mascara=S(i-1:i+1,j-1:j+1);
         
            if esMinLoc(mascara)
                Qz(i-1,j-1) = 1;
                ceros(k,:)=[i,j];
                k=k+1;
                saltProx=1;
            end
            
        end
        
        
    end
end

ceros(:,1)=ceros(:,1)-1;
ceros(:,2)=ceros(:,2)-1;

% ceros(ceros(:,1)>N1,:)=[];
% ceros(ceros(:,2)>M1,:)=[];
% 
% ceros(ceros(:,1)>N1,:)=[];
% ceros(ceros(:,2)>M1,:)=[];

end

function flag = esMinLoc(mascara)
flag=0;
[m,loc]=min(mascara(:));


if loc(1)==5
    flag=1;
end



end