%This function implement the inverse model developed to detect corrodded
%rebar embedded inside concrete or asphalt. To use this function, pass the
%following parameter:
%     Intermediate frequency (IF) signal: E_if
%     Amplitude of the transmitted EM wave: E0
%     Sampling frequency of the system: fs
%     Starting frequency of the FMCW system: F0
%     Bandwidth of the FMCW system: Bc
%     Sweep time of the FMCW system: Tc
function para_estimate = inverse_model(E_if, E0, fs, F0, Bc, Tc)

    n = length(E_if);
    t = (0:(n-1))/fs;

    %Fourier transform to find the frequency of the IF signal
    ft_sig = fft(E_if);
    ps = abs(ft_sig)/n;
    frequency = fs*(0:(n/2-1))/n;

    %Find peaks
    [pks, idx] = findpeaks(ps(1:(n/2)));

    %Find maximum peaks
    [k_peaks_max, idx_peaks_max] = maxk(pks,2);

    chossen_peaks = k_peaks_max;
    chossen_peaks_idx = idx(idx_peaks_max);

    %If the second peak has frequency greater than 7 Hz, make it 7 Hz.
    %This is hard-coded because the 7 Hz peak is covered by the sidelobe of the first peak.
    %This should be improved in the future.
    if chossen_peaks_idx(2) > 7
        chossen_peaks(2) = ps(7);
        chossen_peaks_idx(2) = 7;
    end

    k = length(chossen_peaks);

    peaks_freq = fs*(chossen_peaks_idx-1)/n;

    tau_k = peaks_freq*Tc/Bc;

    S_hat = zeros(n, k);

    for i = 1:k
        temp = cos(2*pi*F0*tau_k(i)+pi*Bc*tau_k(i)/Tc*(2*t-tau_k(i)));
        S_hat(:,i) = temp';
    end

    A_hat = (inv(S_hat'*S_hat)*S_hat'*E_if')';

    A_hat = A_hat./abs(A_hat).*chossen_peaks*2;

    %%Estimate permittivity and layers’ thickness
    simp_fwmodel = @(p) simp_fw_model(p, A_hat, tau_k, E0);
    p0 = [ones(1, k) zeros(1, k)];
    p_hat = fsolve(simp_fwmodel,p0,optimoptions('fsolve','Display','none'));

    epsilon_opt = real([1 p_hat(1:k)]);
    d_opt = real(p_hat(k+1:end));

    %%Applying Levenberg-Marquardt curve fitting algorithm with multiple starting guess
    rng default % For reproducibility

    mu_guess = ones(size(epsilon_opt));
    mu_guess(end) = 1000;

    sigma_guess = zeros(size(epsilon_opt));
    sigma_guess(end) = 1e6;

    para_guess = [mu_guess epsilon_opt sigma_guess d_opt];

    mu_lb = ones(size(mu_guess));
    mu_ub = ones(size(mu_guess))*100;
    mu_ub(end) = 1e5;

    epsilon_lb = ones(size(epsilon_opt));
    epsilon_ub = epsilon_lb*85;
    % epsilon_ub(1) = 1.5;
    epsilon_ub(end) = 1e6;

    sigma_lb = zeros(size(sigma_guess));
    sigma_ub = sigma_lb + 100;
    sigma_ub(end) = 7e7;

    d_lb = d_opt - 0.1;
    d_ub = d_opt + 0.1;

    lb = [mu_lb epsilon_lb sigma_lb d_lb];
    ub = [mu_ub epsilon_ub sigma_ub d_ub];

    fw_func = @(p, t) fw_model(p,t,E0,F0,Bc,Tc);

    problem = createOptimProblem('lsqcurvefit','x0',para_guess,'objective',fw_func,...
        'lb',lb,'ub',ub,'xdata',t,'ydata',E_if);

    ms = MultiStart('Display','none');
    para_estimate = run(ms,problem,50);
    
    plot(t, E_if, "Color","g");
    title("Comparision of real data and generated data from inverse model");
    hold on;
    plot(t, fw_model(para_estimate, t, E0, F0, Bc, Tc), "Color","r", "LineWidth",1);
    xlabel("time (s)");
    ylabel("Amplitude (V/m)");
    legend('Real data','Generated data after estimation')
    hold off;
end

