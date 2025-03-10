function [weighted_left_value, weighted_right_value, f_weighted, p_combined, peak_value, peak_freq, sigma,f_combined_valid,consum,consum_total] = calculate_weighted_iqr(f_values, p_values)
    % f_values 是一个包含所有波段频率数据的 cell 数组
    % p_values 是一个包含所有波段功率数据的 cell 数组

    % 初始化联合功率谱，不进行预设值，直接逐波段相乘
    p_combined = [];  % 空值开始累乘，避免初始值偏差
    f_weighted = f_values{1};  % 假设所有波段的频率范围相同

    % 遍历每个波段数据，使用最大值归一化后进行累乘
    for i = 1:length(p_values)
        f = f_values{i};
        p = p_values{i};

        % 使用最大值归一化
        p_normalized = p / max(p);

        if isempty(p_combined)
            % 第一个波段数据直接赋值
            p_combined = p_normalized;
        else
            % 逐波段累乘
            p_combined = p_combined .* p_normalized;
        end
    end

    % 计算联合功率谱的峰值及其对应的频率
    [peak_value, peak_index] = max(p_combined);
    peak_freq = f_weighted(peak_index);

    % 计算 20% 的噪声阈值
    noise_threshold_combined = peak_value * 0.2;
    
    % 筛选出大于 20% 阈值的数据部分
    valid_indices = p_combined > noise_threshold_combined;
    f_combined_valid = f_weighted(valid_indices);
    p_combined_valid = p_combined(valid_indices);

    % 检查是否有有效数据
    if isempty(p_combined_valid)
        warning('No data above the 20%% noise threshold. Skipping IQR calculation.');
        weighted_left_value = NaN;
        weighted_right_value = NaN;
        sigma = NaN;
        peak_freq = NaN;
        return;
    end

    % 计算 IQR
    consum = cumsum(p_combined_valid);
    consum_total = consum(end);
    percent_25 = quantile(consum, 0.25);
    percent_75 = quantile(consum, 0.75);

    % 找到对应的频率边界值
    left_index = find(consum >= percent_25, 1);
    right_index = find(consum >= percent_75, 1);
    
    weighted_left_value = f_combined_valid(left_index);
    weighted_right_value = f_combined_valid(right_index);
    
    % 计算 sigma 并打印到命令窗口
    sigma = (weighted_right_value - weighted_left_value) / 2;
    fprintf('Peak Height: %.2f,  Sigma: %.2f\n', peak_freq, sigma);

    % 输出最终的联合功率谱结果
    p_combined = p_combined;  % 返回联合后的功率谱



end
