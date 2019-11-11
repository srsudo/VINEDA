
function [CF,y,e,d,fsp]=VINEDA_25oct2018(x,fs,flow,fhigh,Nfb,Dmin,Dmax,Ndb,beta)

% Input:    x: Data samples 
%           fs: Original sampling frequency (Hz)
%           flow: Lower cutoff frequency of the bandpass filter (Hz)
%           fhigh: Upper cutoff frequency of the bandpass filter (Hz)
%           Nfb: Number of frequency bands 
%           Dmin: Minimum duration of the explosions to be detected (s)
%           Dmax: Maximum duration of the explosions to be detected (s)
%           Ndb: Number of duration bands 
%           beta: Factor for reduction of non-stationary noises

% Output:   CF: Characteristic function
%           fsp: Sampling frequency of the characteristic function


% Low pass filter + decimation (the sampling frequency is fsp (Hz) for now on)
x=x-mean(x);
fsp=ceil(2*fhigh);          
xr=resample(x,fsp,fs);   

% Adaptive detrending filter for removing trends (avoiding noises associated with wind, etc)  
B=fir1(round(Dmax*fsp),flow/(fsp/2));      
trend=filtfilt(B,1,xr);
xrd=xr-trend;

% Median filter (avoiding spiky signals, instrumental noise, etc)
y=medfilt1(xrd,round((Dmin/10)*fsp));

% Multiband filter stage + envelope detector (minimizing the effect of
% background noise and providing no edge delay)
th=0:1/fsp:1/flow;
haux=linspace(1,0,length(th));
k=1:Nfb;
fl=flow + (k-1)*(fhigh-flow)/Nfb;
fh=fl + (fhigh-flow)/Nfb;
v=zeros(length(y),Nfb);
ei=zeros(length(y),Nfb);
for k=1:Nfb
    fc=(fl(k)+fh(k))/2;           
    h1=haux.*cos(2*pi*fc*th);   
    v(:,k)=filter(h1,1,y);
    ei(:,k)=abs(hilbert(v(:,k)));
end
e=sum(ei,2);                     

% Discriminant detector (enhancing signals with a sharp rise and gradual
% decay, whilst mitigating the effect of both emerging and impulsive noises
D=linspace(Dmin,Dmax,Ndb);
D=D+1/flow;
di=zeros(length(e),Ndb);
for k=1:Ndb
    LB=round(D(k)*fsp);
    h2=zeros(1,2*LB);
    h2(1:LB)=(1:LB);
    h2(LB+(1:LB))=beta*((1:LB)-(LB+1));
    h2=h2/sqrt(sum(h2.^2));
    di(:,k)=filter(h2,1,e);
    M=floor(length(h2)/2);
    di(:,k)=[di(M+1:end,k) ; zeros(M,1)];
    di(1:M,k)=0;    
end
di(di<=0)=0;
d=prod(di,2);

% Normalization 
V=2^23;     % Normalization value. Maximum value of the input signal x (24-bit digitizers)
mu=1;      % µ-law parameter of the compander
CF=10 * (log10(1+mu*d/V)/log10(1+mu));   % VINEDA law 



