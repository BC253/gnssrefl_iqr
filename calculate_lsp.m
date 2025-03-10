function [f, p] = calculate_lsp(data_matrix, system_band)
    % Calculate Lomb-Scargle Periodogram for given data.
    %
    % Parameters:
    % data_matrix : matrix
    %     Input data where the first column is elevation and the second is SNR.
    % system_band : char
    %     String to indicate the system and band (e.g., 'L1_G', 'E1', 'B1I', etc.).
    %     Supported options: L1_G, L2_G, L5_G, E1, E5a, E5b, E5, E6, B1C, B1I, B3I, B2a, B2I, B2b, B2a+b, L1_R, L2_R.
    %
    % Returns:
    % f : vector
    %     Frequency values.
    % p : vector
    %     Power values of the Lomb-Scargle periodogram.
    % Define frequency values based on system_band

    frequency_dict = struct('L1_G', 1575.42, ... % GPS L1
                            'L2_G', 1227.60, ... % GPS L2
                            'L5_G', 1176.45, ... % GPS L5
                            'L1_E', 1575.42, ...   % Galileo E1
                            'L5_E', 1176.45, ...  % Galileo E5a
                            'L7_E', 1207.14, ...  % Galileo E5b
                            'L8_E', 1191.795, ...  % Galileo E5 (E5a + E5b)
                            'L6_E', 1278.75, ...   % Galileo E6
                            'L2_C', 1561.098, ...  % BDS B1C
                            'L7_C', 1207.14, ... % BDS B1I
                            'L6_C', 1268.52, ...  % BDS B3I
                            'L1_C', 1575.42, ...  % BDS B2a
                            'L5_C', 1176.45, ...  % BDS B2I
                            'L9_C', 1207.14, ...  % BDS B2b
                            'L8_C', 1191.795, ... % BDS B2a+b
                            'L1_R', 1600, ...    % GLONASS L1 (Approximate)
                            'L2_R', 1250);       % GLONASS L2 (Approximate)

    if ~isfield(frequency_dict, system_band)
        error('Unsupported system_band: %s. Please use one of %s', system_band, strjoin(fieldnames(frequency_dict), ', '));
    end

    % Get frequency based on the selected system band
    frequency_mhz = frequency_dict.(system_band);

    % Calculate constants
    c = 299792458;  % Speed of light (m/s)
    cf = c / (frequency_mhz * 1e6) / 2;

    % Extract data from the matrix
    e = data_matrix(:, 1);  % Elevation angles
    snr = data_matrix(:, 2);  % SNR values

    % Calculate Lomb-Scargle Periodogram parameters
    minRH = 0;
    maxRH = 6;
    desiredPrecision = 0.01;
    [ofac, hifac] = get_ofac_hifac(e, cf, maxRH, desiredPrecision);

    % Prepare data for Lomb-Scargle
    sineE = sind(e);
    [sortedX, j] = sort(sineE);
    sortedY = snr(j);

    % Calculate Lomb-Scargle Periodogram
    [f, p] = lomb(sortedX / cf, sortedY, ofac, hifac);

    % Filter frequencies greater than minRH
    valid_indices = f > minRH;
    f = f(valid_indices);
    p = p(valid_indices);
    if ~isfield(frequency_dict, system_band)
    warning('Unsupported system_band: %s. Skipping calculation for this band.', system_band);
    f = [];
    p = [];
    return;
    
end

% Note: Functions like get_ofac_hifac and lomb are assumed to be implemented elsewhere.