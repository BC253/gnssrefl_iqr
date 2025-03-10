function calculate_and_plot_iqr(f, p, sat, system_band_combined)
    % 找到最高峰值并计算20%噪声阈值 (Find peak value and calculate 20% noise threshold)
    [peak_value, peak_index] = max(p);
    noise_threshold = peak_value * 0.2;

    % 打印峰值 (Print peak value)
    fprintf('Peak Value: %.2f\n', f(peak_index));


    % 去除噪声部分 (Remove noise portion)
    valid_indices = p > noise_threshold;
    f_valid = f(valid_indices);
    p_valid = p(valid_indices);

    % 初始化峰值范围变量
    f_peak_range = [];
    p_peak_range = [];

    % 检查峰值数量 (Check number of peaks)
    [peaks, peak_locs] = findpeaks(p); % 找出所有峰值
    main_peak_index = find(peak_locs == peak_index);

    % 找到主峰左边的第一个峰值
    left_peak_index = main_peak_index - 1;
    if left_peak_index > 0
        left_peak = peak_locs(left_peak_index);
    else
        left_peak = []; % 没有左峰值
    end

    % 找到主峰右边的第一个峰值
    right_peak_index = main_peak_index + 1;
    if right_peak_index <= length(peaks)
        right_peak = peak_locs(right_peak_index);
    else
        right_peak = []; % 没有右峰值
    end

    % 根据左右峰值的存在情况确定范围
    if ~isempty(left_peak) && ~isempty(right_peak)
        % 左右都有峰值，使用左右峰值之间的范围
        f_peak_range = f(left_peak:right_peak);
        p_peak_range = p(left_peak:right_peak);
    elseif ~isempty(left_peak)
        % 只有左峰值，使用左峰值到主峰右侧20%以上的范围
        f_peak_range = f(left_peak:peak_index);
        p_peak_range = p(left_peak:peak_index);
    elseif ~isempty(right_peak)
        % 只有右峰值，使用主峰左侧20%以上到右峰值的范围
        f_peak_range = f(peak_index:right_peak);
        p_peak_range = p(peak_index:right_peak);
    else
        % 没有左右峰值，只使用20%以上的范围
        f_peak_range = f_valid;
        p_peak_range = p_valid;
    end

    % 计算累加范围内的值 (Cumulative sum within range)
    consum = cumsum(p_peak_range);
    consum_total = consum(end);

    percent_25 = quantile(consum, 0.25);
    percent_75 = quantile(consum, 0.75);

    % 找到对应的频率值 (Find corresponding frequency values)
    left_index = find(consum >= percent_25, 1);
    right_index = find(consum >= percent_75, 1);
        % 计算sigma值 (Calculate sigma value)
    weighted_left_value = f_peak_range(left_index);
    weighted_right_value = f_peak_range(right_index);
    sigma_value = (weighted_right_value - weighted_left_value) / 2;
    fprintf('Sigma Value: %.2f', sigma_value);

    % 确保找到边界值 (Ensure boundary values are found)
    if isempty(left_index) || isempty(right_index)
        warning('无法找到25%或75%的频率值');
        plot(f, p, 'DisplayName', ['Satellite ' num2str(sat) ' - ' system_band_combined], 'Color', [0.7, 0.7, 0.7]);
        xlabel('Frequency');
        ylabel('Power');
        title(['Satellite ' num2str(sat) ' - ' system_band_combined]);
        grid on;
        return;
    end

    left_value = f_peak_range(left_index);
    right_value = f_peak_range(right_index);

    % 将数据分为三个部分: 低于0.25范围、中间0.25-0.75范围、高于0.75范围
    lower_range_indices = f < left_value;
    middle_range_indices = f >= left_value & f <= right_value;
    upper_range_indices = f > right_value;

    % 绘制不同范围的数据
    plot(f(lower_range_indices), p(lower_range_indices), 'Color', [0.7, 0.7, 0.7],'HandleVisibility', 'off');
    hold on;
    plot(f(middle_range_indices), p(middle_range_indices), 'Color', [0.1, 0.6, 0.1], 'LineWidth', 1.5, 'DisplayName', '25%-75% Range');
    plot(f(upper_range_indices), p(upper_range_indices), 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'original data');

    % 添加噪声阈值线 (Add noise threshold line)
    yline(noise_threshold, '--r', ['20% Noise: ' num2str(noise_threshold, '%.2f')], 'LabelVerticalAlignment', 'bottom','DisplayName', '20% Noise');

    % 绘制IQR边界线 (Draw IQR boundary lines)
    xline(left_value, '--b', ['25% Value: ' num2str(left_value, '%.2f')], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left','DisplayName', '25% Value');
    xline(right_value, '--r', ['75% Value: ' num2str(right_value, '%.2f')], 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right','DisplayName', '75% Value');

    % 添加峰值线 (Add peak value line)
    xline(f(peak_index), '--k', ['Peak Value: ' num2str(f(peak_index), '%.2f')], 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right','DisplayName', 'Peak Value');

    % 设置图例和标题 (Set legend and title)
    xlabel('Reflection height [m]');
    ylabel('Voltage');
    title(['Satellite ' num2str(sat) ' - ' system_band_combined]);
    legend;
    grid on;
    hold off;
     figure 
        y = consum/consum_total;
        plot(f_peak_range,y)
        hold on
        index1 = find(f_peak_range == left_value);
        index2 = find(f_peak_range == right_value);
        yline(y(index1), '--b', ['25% Value: ' num2str(left_value, '%.2f')], ...
               'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'DisplayName', '25% Boundary');
        
        yline(y(index2), '--r', ['75% Value: ' num2str(right_value, '%.2f')], ...
               'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'DisplayName', '75% Boundary');
        title('IQR');
        xlabel('Reflector Height [m]');
        ylabel('Normalized Power');
        grid on;
        hold off
end
