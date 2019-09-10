
function [CF,fsp]=VINEDA(x,fs,flow,fhigh,L,Nc,beta)

% Input:    x: Data samples 
%           fs: Original sampling frequency (Hz)
%           flow: Lower cutoff frequency of the bandpass filter (Hz)
%           fhigh: Upper cutoff frequency of the bandpass filter (Hz)
%           L: Average duration of the explosions to be detected (s)
%           Nc: Number of subbands used in the filter bank
%           beta: Penalty factor por non-impulsive onsets

% Output:   CF: Characteristic function
%           fsp: Sampling frequency of the characteristic function


% Low pass filter + decimation (the sampling frequency is fsp (Hz) for now on)
x=x-mean(x);
fsp=ceil(2*fhigh); 
xr=resample(x,fsp,fs);   

% Adaptive detrending filter for removing trends (avoiding noises associated with wind, etc)  
B=fir1(1*L*fsp,flow/(fsp/2));            
trend=filtfilt(B,1,xr);
xrd=xr-trend;

% Median filter (avoiding spiky signals, instrumental noise, etc)
y=medfilt1(xrd,round((L/10)*fsp));

% Multi-band filter stage + envelope detector (minimizing the effect of
% background noise and providing no edge delay)
lh=round(L*fsp);
th=0:1/fsp:(lh-1)/fsp;
haux=linspace(1,0,lh);
k=1:Nc;
fl=flow + (k-1)*(fhigh-flow)/Nc;
fh=fl + (fhigh-flow)/Nc;
v=zeros(length(y),Nc);
ei=zeros(length(y),Nc);
for k=1:Nc,
    fc=(fl(k)+fh(k))/2;           
    h1=haux.*sin(2*pi*fc*th);    
    v(:,k)=filter(h1,1,y);
    ei(:,k)=abs(hilbert(v(:,k))); 
end
e=sum(ei,2);                      

% Discriminant detector (enhancing signals with a sharp rise and gradual
% decay with a duration in the order of L secs, whilst mitigating the effect of non-stationary noises)
LB=round(L*fsp);
h2(1:LB)=(1:LB); 
h2(LB+(1:LB))=beta*((1:LB)-(LB+1)); 
h2=h2/sqrt(sum(h2.^2));
CF=filter(h2,1,e);
M=floor(length(h2)/2);             
CF=[CF(M+1:end) ; zeros(M,1)];     
CF(1:M)=0;
CF(CF<=0)=0;

% Normalization using a 30-secs window for estimating background noise level
N=round(30*fsp);            
h3=(1/N)*ones(1,N);
n=filter(h3,1,e);
M=floor(N/2);
n=[n(M+1:end) ; zeros(M,1)];      
n(1:M)=n(M+1);
n(end-M+1:end)=n(end-M);
CF(CF~=0)=CF(CF~=0)./n(CF~=0);         

