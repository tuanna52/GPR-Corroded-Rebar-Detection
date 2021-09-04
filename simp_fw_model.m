%This function implement the simplified forward model developed used in
%inverse model to estimate optimized values for layers' permittivity and thickness.
%To use this function, pass the following parameter:
%     Input parameters: [permittivity_array  layers_thickness_array]
%     Estimated amplitude of k elements of the IF signal: A-hat
%     Two-way time travel from kth interface to the antenna: tau_k
%     Amplitude of the transmitted EM wave: E0
function output = simp_fw_model(input_para, A_hat, tau_k, E0)
    c0 = 299792458; %Velocity of light in vacuum
    
    para_num = length(A_hat);
    output = zeros(1, para_num*2);

    epsilon_hat = [1 input_para(1:para_num)];
    d_hat = input_para(para_num+1:end);

    R_coef_fw = ones(size(d_hat)); %Reflection coefficient from k-th layer to (k+1)-th layer
    R_coef_bw = ones(size(d_hat)); %Reflection coefficient from (k+1)-th layer to k-th layer
    T_coef_fw = ones(size(d_hat)); %Transmission coefficient from k-th layer to (k+1)-th layer
    T_coef_bw = ones(size(d_hat)); %Transmission coefficient from (k+1)-th layer to k-th layer

    for k = 1:para_num
        R_coef_fw(k) = (sqrt(epsilon_hat(k+1))-sqrt(epsilon_hat(k)))/(sqrt(epsilon_hat(k+1))+sqrt(epsilon_hat(k)));
        R_coef_bw(k) = (sqrt(epsilon_hat(k))-sqrt(epsilon_hat(k+1)))/(sqrt(epsilon_hat(k+1))+sqrt(epsilon_hat(k)));
        T_coef_fw(k) = (2*sqrt(epsilon_hat(k+1)))/(sqrt(epsilon_hat(k+1))+sqrt(epsilon_hat(k)));
        T_coef_bw(k) = (2*sqrt(epsilon_hat(k)))/(sqrt(epsilon_hat(k+1))+sqrt(epsilon_hat(k)));
    end

    tz = 0;
    for k = 1:para_num
        tz = tz + 2*d_hat(k)*sqrt(epsilon_hat(k))/c0;
        output(k) = tz - tau_k(k);
    end

    output(1+para_num) = E0^2/2*R_coef_fw(1) - A_hat(1);

    T_coef = 1;
    for k = 2:para_num
        T_coef = T_coef*T_coef_fw(k-1)*T_coef_bw(k-1);
        output(k+para_num) = E0^2/2*R_coef_fw(k)*T_coef - A_hat(k);
    end
    
end

