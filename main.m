% Wireless Communication Systems_Final Project
% Due date: 20190119
% 107064522

clear all
close all
%% Input parameters
N = 4;  % The number of antennas
delta = 0.5;    % The normalized antenna separation

% for comparison
% % % N = 16;
% % % delta = 1/2;
%% Output results
%% The angular-domain radiation/reception basis
L = delta*N; % normalized length of the antenna array
U = zeros(N, N);
k = 1:N;
for i = 1:N
    Omega = (k-1)/L;
    U(:,i) = exp(-j*2*pi*(i-1)*delta*Omega)/sqrt(N);    % (p.206)
    %fprintf('%f, %f\n',i-1,Omega);
end
U;

%% The correlation between different basis vectors
COR = zeros(N, N);
for i = 1:N
    for k = 1:N
        COR(i,k) = U(:,i)'*U(:,k);  % correlation of vectors
    end
end
COR = abs(COR);
COR;

% comparison
% % % N = 16;
% % % delta = 1/8;

Omega = -2:0.001:2; 
corr = abs((1/N)*exp(j*pi*Omega*delta*(N-1)).*sin(pi*N*delta*Omega)./sin(pi*delta*Omega)); %(7.59)
plot(Omega,corr);
%% The gain pattern of the ULA
phi_0 = pi/2; % the angle of incidence of the LoS onto the Rx antenna array
omega = cos(phi_0);   % directional cosine
phi = 0:0.01:2*pi;

% comparison
% % % phi_0 = 0;
% % % N = 16;
% % % delta = 1/2;
gain = zeros(N,length(phi));
for i = 1:N
    Omega = cos(phi)-cos(phi_0);    % (7.42)
    %Omega = cos(phi)-i/L;    % (7.42)
    gain(i,:) = abs((1/N)*exp(j*pi*Omega*delta*(N-1)).*sin(pi*N*delta*Omega)./sin(pi*delta*Omega)); %(7.59)
    %figure
    %polar(phi, gain(i,:));
end

figure;
polar(phi, gain(1,:))
hold on
for i = 2:N
    polar(phi, gain(i,:))
end
hold off

%% The gain of the desired signal for using different radiation/reception beams
x = 1;  % attenuation of the path
d = 10^6;   % distance between Tx and first Rx antennas
lambda = (3*10^8)/(10^9);   % speed of light/carrier frequency

d_i = zeros(1,N);
for i = 1:N
    d_i(i) = d+(i-1)*delta*lambda*omega;	% distance between antennas (7.19)
end
h = x*exp(-j*2*pi*d_i/lambda); % The channel gain of desired signal (7.17)
g = zeros(1,N);
for i = 1:N
    g(i) = sum(U(:,i).*h');
end
g = abs(g);
g;

%% The signal-to-interference power ratio (SIR) for using different beams
% for interference condition, we update d and h
phi_interference = pi/2;
d_interference = zeros(1,N);    % distance between antennas (interference version)
for i = 1:N
    d_interference(i) = d+(i-1)*delta*lambda*cos(phi_interference);	% distance between antennas (7.19)
end

h_interference = x*exp(-j*2*pi*d_interference/lambda);
g_interference = zeros(1,N);
for i = 1:N
    g_interference(i) = sum(U(:,i).*h_interference');
end
g_interference = abs(g_interference);

% computing SIR (signal power/interference power)
P_signal = g.^2;
P_interference = g_interference.^2;
SIR = P_signal./P_interference;
SIR_dB = 10*log10(SIR);

%% The SIR of multiple input signals (multiple reception directions) with diversity combining (considering fading for the signals and interference)
% generating QPSK signals
no_samples = 10^5;   % Number of bits to be transmitted
E = 4;
s = round(rand(size(zeros(1,no_samples))));
symbol = zeros(1,no_samples/2);
for i = 1:no_samples/2
    if(~s(2*i))
        x = sqrt(E/2);
    else
        x = -sqrt(E/2);
    end
    if(~s(2*i-1))
        y = sqrt(E/2);
    else
        y = -sqrt(E/2);
    end
    symbol(i) = x+j*y;
end
symbol = repmat(symbol, N, 1);

% through Rayleigh fading channel
I_data = normrnd(0,1^0.5,size(symbol));   % Inphase bipolar seq.
Q_data = normrnd(0,1^0.5,size(symbol));   % Quadrature bipolar seq.
channel = I_data+j*Q_data;   % QPSK signal
channel_i = channel;

% SC
branch = zeros(size(symbol));
for i = 1:length(h)
    branch(i,:) = symbol(i,:).*conj(channel(i,:)*h(i));   % phase elimination
end
tmp1 = zeros(size(symbol));
tmp2 = tmp1;
for i = 1:length(h)                      
    tmp1(i,:) = h(i)*channel(i,:).* branch(i,:);
    tmp2(i,:) = h_interference(i)*channel(i,:).* branch(i,:);
end
SC = max(sum(tmp1,1));
SC_interference = max(sum(tmp2,1));

SC_SIR = sum(abs(SC(1,:)).^2)/sum(abs(SC_interference(1,:)).^2); %signal power/interference power
SC_SIR_dB = 10*log10(SC_SIR);	%signal power/interference power (dB scale)


% MRC
branch = zeros(size(symbol));
for i = 1:length(h)
    branch(i,:) = symbol(i,:).*conj(channel(i,:)*h(i));   % phase elimination
end
tmp1 = zeros(size(symbol));
tmp2 = tmp1;
for i = 1:length(h)                      
    tmp1(i,:) = h(i)*channel(i,:).* branch(i,:);
    tmp2(i,:) = h_interference(i)*channel(i,:).* branch(i,:);
end
MRC = sum(tmp1,1);
MRC_interference = sum(tmp2,1);

MRC_SIR = sum(abs(MRC(1,:)).^2)/sum(abs(MRC_interference(1,:)).^2); %signal power/interference power
MRC_SIR_dB = 10*log10(MRC_SIR);	%signal power/interference power (dB scale)

% EGC
branch = zeros(size(symbol));
for i = 1:length(h)
    branch(i,:) = symbol(i,:).*conj(channel(i,:)*h(i))./abs(channel(i,:)*h(i));   % phase elimination   
end

tmp1 = zeros(size(symbol));
tmp2 = tmp1;
for i = 1:length(h)                      
    tmp1(i,:) = h(i)*channel(i,:).* branch(i,:);
    tmp2(i,:) = h_interference(i)*channel(i,:).* branch(1,:); 
end
EGC = sum(tmp1,1);
EGC_interference = sum(tmp2,1);

EGC_SIR = sum(abs(EGC(1,:)).^2)/sum(abs(EGC_interference(1,:)).^2); %signal power/interference power
EGC_SIR_dB = 10*log10(EGC_SIR);	%signal power/interference power (dB scale)