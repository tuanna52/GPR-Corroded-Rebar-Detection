%This function implement the forward model developed to detect corrodded
%rebar embedded inside concrete or asphalt. To use this function, pass the
%following parameter:
%     Input parameters: [permeability_array permittivity_array conductivity_array layers_thickness_array]
%     Recorded time: Time array from the beginning to the end
%     Amplitude of the transmitted EM wave: E0
%     Starting frequency of the FMCW system: F0
%     Bandwidth of the FMCW system: Bc
%     Sweep time of the FMCW system: Tc
function [output] = fw_model(input_param, t, E0, F0, Bc, Tc)
%%  Input  
    param_number = length(input_param);

    mu0 = 4*pi*1e-7; %Vacuum permeability
    mu = input_param(1:((param_number+1)/4)) * mu0; %Permeability of medium layers

    epsilon0 = 8.8541853e-12; %Vacuum permittivity
    epsilon = input_param(((param_number+1)/4+1):((param_number+1)/4*2)) * epsilon0; %Permittivity of medium layers

    sigma = input_param((((param_number+1)/4)*2+1):(((param_number+1)/4)*3)); %Conductivity of medium layers

    d = input_param((((param_number+1)/4)*3+1):param_number); %Layers' thichness (m)
%%  Derived parameters

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
    
%%  Forward model

    Ek = E0*R_coef_fw(1);
    tau = 2*d(1)/vp(1);
    E_if = E0*abs(Ek)/2*cos(2*pi*F0*tau + pi*Bc*tau/Tc*(2*t-tau) - angle(Ek));

    Ek = 1;
    T_coef = 1;
    tz = 0;
    for k = 2:length(d)
        tz = tz + 2*d(k-1)/vp(k-1);
        T_coef = T_coef*T_coef_fw(k-1)*T_coef_bw(k-1);
        for q = 1:2 %Assume there is two multiple reflections inside each medium
            tau = 2*q*d(k)/vp(k) + tz;
            Ek = E0*(R_coef_fw(k)^q)*(R_coef_bw(k-1)^(q-1))*T_coef;
            E_temp = E0*abs(Ek)/2*cos(2*pi*F0*tau + pi*Bc*tau/Tc*(2*t-tau) - angle(Ek));
            E_if = E_if + E_temp;
        end
    end

    output = E_if;
end

