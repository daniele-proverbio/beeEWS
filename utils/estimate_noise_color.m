% Color of Noise Estimation via PSD Slope Fitting
% ================================================
% Fits 1/f^beta to estimate noise color:
%   beta = 0 → white noise
%   beta = 1 → pink (1/f) noise
%   beta = 2 → red/brown noise
%   beta = -1 → blue noise
%   beta = -2 → violet noise

function noise_color = estimate_noise_color(signal, fs)
    % INPUTS: signal (vector of time series) and fs (sampling frequency in
    % Hz (default: 1))
    % OUTPUT: noise_color - struct with beta, R², color label, and PSD data

    if nargin < 2
        fs = 1;
    end

    N = length(signal);

    % Compute Power Spectral Density via Welch's method
    %     Use ~8 segments with 50% overlap for stable estimate
    seg_len   = floor(N / 4.5);          % segment length
    n_overlap = floor(seg_len * 0.5);    % 50% overlap
    n_fft     = 2^nextpow2(seg_len);     % FFT length (power of 2)

    [psd, freq] = pwelch(signal, hann(seg_len), n_overlap, n_fft, fs);

    % Remove DC component and Nyquist (boundary artifacts)
    freq(1) = [];
    psd(1)  = [];

    % Convert to log-log scale for linear fitting
    log_freq = log10(freq);
    log_psd  = log10(psd);

    % Fit a line in log-log space:  log(PSD) = -beta * log(f) + const
    %     Use only the central 80% of frequencies to avoid edge effects
    n        = length(log_freq);
    idx_fit  = round(0.1*n) : round(0.9*n);

    p        = polyfit(log_freq(idx_fit), log_psd(idx_fit), 1);
    beta     = -p(1);       % slope of PSD ∝ 1/f^beta  →  beta = -slope
    offset   =  p(2);

    % Compute R2 (goodness of fit)
    psd_fit      = polyval(p, log_freq(idx_fit));
    ss_res       = sum((log_psd(idx_fit) - psd_fit).^2);
    ss_tot       = sum((log_psd(idx_fit) - mean(log_psd(idx_fit))).^2);
    r_squared    = 1 - ss_res / ss_tot;

    % Classify noise color
    color_label = classify_noise_color(beta);

    % Package results
    noise_color.beta        = beta;
    noise_color.r_squared   = r_squared;
    noise_color.color       = color_label;
    noise_color.freq        = freq;
    noise_color.psd         = psd;
    noise_color.fit_slope   = p(1);
    noise_color.fit_offset  = offset;
    noise_color.idx_fit     = idx_fit;

    % % Print summary
    % fprintf('=== Noise Color Estimation ===\n');
    % fprintf('  Spectral exponent (beta) : %.3f\n', beta);
    % fprintf('  Goodness of fit (R²)     : %.4f\n', r_squared);
    % fprintf('  Noise color              : %s\n\n', color_label);

    
end


%% Helper functions 

% Label noise color
function label = classify_noise_color(beta)
    if     beta < -1.5,          label = 'Violet (β ≈ -2)';
    elseif beta < -0.5,          label = 'Blue (β ≈ -1)';
    elseif abs(beta) <= 0.5,     label = 'White (β ≈ 0)';
    elseif beta <= 1.5,          label = 'Pink / 1/f (β ≈ 1)';
    elseif beta <= 2.5,          label = 'Red / Brownian (β ≈ 2)';
    else,                        label = 'Black (β > 2)';
    end
end


