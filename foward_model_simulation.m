%% Input parameters
E0 = 1; %Amplitude of transmitted signal

mu0 = 4*pi*1e-7; %Vacuum permeability
mu = [1 1 100] * mu0; %Permeability of medium layers

epsilon0 = 8.8541853e-12; %Vacuum permittivity
epsilon = [1 6 1000] * epsilon0; %Permittivity of medium layers

sigma = [1e-9 5e-3 7e6]; %Conductivity of medium layers

d = [10e-2 15e-2]; %Layers' thichness (m)

F0 = 24e9; %Start frequency of the chirp (Hz)
Bc = 1.5e9; %Bandwidth of the chirp (Hz)
Tc = 300e-3; %Sweep time of the chirp (s)

SNR = 10; %Signal to noise ratio

fs = 5000; %Sampling frequency of the system (Hz)
obs_time = 3*Tc; %Observation time (s)
t = 0:1/fs:(obs_time - 1/fs);

%% Derived parameters
omega_c = 2*pi*(F0+Bc/2); %Angle velocity of center frequency

gamma = ones(size(mu)); %Propagation constant of medium
eta = ones(size(mu)); %Intrinsic impedance of medium
vp = ones(size(mu)); %EM Wave velocity in medium

for k = 1:length(mu)
    gamma(k) = sqrt(1i*omega_c*mu(k)*(sigma(k)+1i*omega_c*epsilon(k)));
    eta(k) = sqrt(1i*omega_c*mu(k)/(sigma(k)+1i*omega_c*epsilon(k)));
    vp(k) = 1/sqrt(mu(k)*epsilon(k)/2*(sqrt(1+(sigma(k)/(omega_c*epsilon(k)))^2)+1));
end

R_coef_fw = ones(size(d)); %Reflection coefficient from k-th layer to (k+1)-th layer
R_coef_bw = ones(size(d)); %Reflection coefficient from (k+1)-th layer to k-th layer
T_coef_fw = ones(size(d)); %Transmission coefficient from k-th layer to (k+1)-th layer
T_coef_bw = ones(size(d)); %Transmission coefficient from (k+1)-th layer to k-th layer

z = 0;
for k = 1:length(d)
    z = z + d(k);
    R_coef_fw(k) = exp(-2*gamma(k)*z)*((eta(k+1)-eta(k))/(eta(k+1)+eta(k)));
    R_coef_bw(k) = exp(2*gamma(k+1)*z)*((eta(k)-eta(k+1))/(eta(k+1)+eta(k)));
    T_coef_fw(k) = exp((gamma(k+1)-gamma(k))*z)*2*eta(k+1)/(eta(k+1)+eta(k));
    T_coef_bw(k) = exp((gamma(k+1)-gamma(k))*z)*2*eta(k)/(eta(k+1)+eta(k));
end

%% Forward model
Ek = E0*R_coef_fw(1);
tau = 2*d(1)/vp(1);
E_if = E0*abs(Ek)/2*cos(2*pi*F0*tau + pi*Bc*tau/Tc*(2*t-tau) - angle(Ek));

Ek = 1;
T_coef = 1;
tz = 0;
for k = 2:length(d)
    tz = tz + 2*d(k-1)/vp(k-1);
    T_coef = T_coef*T_coef_fw(k-1)*T_coef_bw(k-1);
    for q = 1:2 %Assume there is two multiple reflections inside each layer
        tau = 2*q*d(k)/vp(k) + tz;
        Ek = E0*(R_coef_fw(k)^q)*(R_coef_bw(k-1)^(q-1))*T_coef;
        E_temp = E0*abs(Ek)/2*cos(2*pi*F0*tau + pi*Bc*tau/Tc*(2*t-tau) - angle(Ek));
        E_if = E_if + E_temp;
    end
end

clear E_temp;

%Final IF signal with white Gaussian noise added
E_if = awgn(E_if, SNR);

plot(t, E_if);

%Fourier transform to find the frequency of the IF signal
ps = abs(fft(E_if))/length(t);
frequency = fs*(0:(length(t)/45))/length(t);

figure;
plot(frequency,ps(1:(length(t)/45)+1));