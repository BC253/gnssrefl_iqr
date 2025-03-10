function plot_gnss_data(filepaths, satellite_numbers, bands, plot_type, save_folder,iqr, error_bar)
    % This function reads data, calculates the LSP for selected satellite(s) and band(s), 
    % and plots the results based on the given plot type.
    %
    % Inputs:
    % - filepaths: cell array of paths to the .txt files containing data.
    % - satellite_numbers: array of satellites to analyze.
    % - bands: cell array of bands of interest (e.g. 'L1', 'L2').
    % - plot_type: type of plot to create ('single_satellite_band', 'same_satellite_different_bands', 
    %   'different_satellites_same_band', 'all_data').
    % - save_folder: folder to save the plot, empty if saving is not needed. 
    % - iqr, true for iqr ,false for normal plot

    % Read and process data from each file in the filepaths
    
    % Extract satellite data and compute LSP
    peak_sigma_data = struct();% This is for Error bar
    satellite_data = struct();
    for i = 1:length(filepaths)
        data_matrix = readmatrix(filepaths{i});
        satellite_num = extract_satellite_number(filepaths{i}); % A function to get satellite number from filename
        band = extract_band(filepaths{i}); % A function to extract the band information from filename
        system = extract_system(filepaths{i});
    
        % Create a new variable for the combined system and band
        system_band_combined = [band, '_', system];
        
        % Store data based on satellite and band
        if ismember(satellite_num, satellite_numbers) && ismember(band, bands)
            if ~isfield(satellite_data, ['sat_' num2str(satellite_num)])
                satellite_data.(['sat_' num2str(satellite_num)]) = struct();
            end
            satellite_data.(['sat_' num2str(satellite_num)]).(system_band_combined) = data_matrix;  
            
        end
    end
    
    % Calculate LSP for each selected satellite and band
    lsp_results = struct();
    for sat = satellite_numbers
        sat_field = ['sat_' num2str(sat)];
        if isfield(satellite_data, sat_field)
            bands_in_data = fieldnames(satellite_data.(sat_field));
            for j = 1:length(bands_in_data)
                band = bands_in_data{j};
                system_band_combined = bands_in_data{j};
                data_matrix = satellite_data.(sat_field).(system_band_combined);
                % Pass the combined system_band to calculate_lsp
                [f, p] = calculate_lsp(data_matrix, system_band_combined);
                lsp_results.(sat_field).(band) = struct('f', f, 'p', p);
                 %  if we use error bar ,calculate peak and sigma
            if error_bar
                    [~, ~, ~, ~,~ , peak_value, sigma, ~, ~, ~] = calculate_weighted_iqr({f}, {p});
                    if ~isfield(peak_sigma_data, system_band_combined)
                        peak_sigma_data.(system_band_combined) = struct('satellite_ids', [], 'peak_vals', [], 'sigma_vals', []);
                    end
                    % 将结果存储在结构数组中
                    peak_sigma_data.(system_band_combined).satellite_ids(end+1) = sat; % 记录卫星 ID
                    peak_sigma_data.(system_band_combined).peak_vals(end+1) = peak_value; % 记录 Peak Value
                    peak_sigma_data.(system_band_combined).sigma_vals(end+1) = sigma; % 记录 Sigma
            
            end
            end
        end
    end




    % Plotting the results based on the specified plot type
    figure_handle = figure;
    hold on;
    switch plot_type
        case 'single_satellite_band'
            % 单卫星单波段绘图 (Single satellite, single band plot)
            plot_title = 'single satellite band';
            if ~isempty(satellite_numbers) && ~isempty(bands)
                sat = satellite_numbers(1);
                band = bands{1};
                sat_field = ['sat_' num2str(sat)];
                if isfield(lsp_results, sat_field) && isfield(lsp_results.(sat_field), system_band_combined)
                    f = lsp_results.(sat_field).(system_band_combined).f;
                    p = lsp_results.(sat_field).(system_band_combined).p;
                    
                    if iqr % 是否计算和绘制 IQR (Calculate and plot IQR if requested)
                       calculate_and_plot_iqr(f, p, sat, system_band_combined);
                    else    
                        % 不受用iqr绘制原始数据 (Plot original data without iqr)
                        plot(f, p, 'DisplayName', ['Satellite ' num2str(sat) ' - ' single_satellite_band]);
                    end

                    % 绘制图例 (Show legend)
                    legend('show');
                end
            end
case 'same_satellite_different_bands'
    % 单卫星多波段绘图 (Single satellite, multiple bands plot)
    plot_title = 'same satellite different bands';
    sat = satellite_numbers(1);
    sat_field = ['sat_' num2str(sat)];
    f_values = {};  % 存储每个波段的频率数据
    p_values = {};  % 存储每个波段的功率数据

    if isfield(lsp_results, sat_field)
        bands_in_data = fieldnames(lsp_results.(sat_field));

        % 遍历每个波段数据
        for j = 1:length(bands_in_data)
            band = bands_in_data{j};
            f = lsp_results.(sat_field).(band).f;
            p = lsp_results.(sat_field).(band).p;

            % 将每个波段的频率和功率数据分别存储
            f_values{end+1} = f;
            p_values{end+1} = p;

            if iqr
                % 绘制灰色版本，用于 IQR 分析
                plot(f, p, 'Color', [0.7, 0.7, 0.7], 'DisplayName', ['Satellite ' num2str(sat) ' - ' band]);
                hold on;

                % 绘制每个波段的20%阈值线
                noise_threshold = max(p) * 0.2;
                yline(noise_threshold, '--r', ['20% Threshold - ' band], 'HandleVisibility', 'off');
            else
                % 如果未请求 IQR，绘制彩色的原始数据
                plot(f, p, 'DisplayName', ['Satellite ' num2str(sat) ' - ' band]);
                hold on;
            end
        end

        % 如果请求 IQR，调用 IQR 计算函数，获取整体的 25% 和 75% 边界范围
        if iqr
            [combined_left_value, combined_right_value, f_weighted, p_combined, peak_value, peak_freq, sigma,f_combined_valid,consum,consum_total] = calculate_weighted_iqr(f_values, p_values);

            % 在第一个图中，将整体的 25%-75% 区间应用到所有波段
            for j = 1:length(bands_in_data)
                band = bands_in_data{j};
                f = f_values{j};
                p = p_values{j};

                % 在整体的 25%-75% 区间内绘制每个波段的数据
                valid_indices = (f >= combined_left_value & f <= combined_right_value);
                plot(f(valid_indices), p(valid_indices), 'Color', [0.2, 0.6, 0.8], 'LineWidth', 1.5, 'DisplayName', ['25%-75% Range - ' band]);
            end

            % 在第一个图中标注整体 25% 和 75% 的边界
            xline(combined_left_value, '--b', ['25% Value: ' num2str(combined_left_value, '%.2f')], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right');
            xline(combined_right_value, '--g', ['75% Value: ' num2str(combined_right_value, '%.2f')], 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
            xline(peak_freq, '--k', ['Peak Frequency: ' num2str(peak_freq, '%.2f')], 'LabelVerticalAlignment', 'top');
            
          
        end
    end

    % 设置图例和标题
    legend('show');
    title(plot_title);
    xlabel('Reflector Height [m]');
    ylabel('Power');
    grid on;











case 'different_satellites_same_band'
    % 不同卫星相同波段绘图 (Different satellites, same band plot)
    plot_title = 'Different Satellites Same Band';
    band = bands{1}; % Select a single band for plotting
    f_values = {};   % 存储每个卫星的频率数据
    p_values = {};   % 存储每个卫星的功率数据

    % 遍历每个卫星的数据 (Iterate over each satellite's data)
    for i = 1:length(satellite_numbers)
        sat = satellite_numbers(i);
        sat_field = ['sat_' num2str(sat)];

        % 构造合适的 system_band_combined
        system_band_combined = [band, '_', extract_system(filepaths{1})];

        if isfield(lsp_results, sat_field) && isfield(lsp_results.(sat_field), system_band_combined)
            % 获取频率和功率数据
            f = lsp_results.(sat_field).(system_band_combined).f;
            p = lsp_results.(sat_field).(system_band_combined).p;

            % 将当前卫星的数据添加到 f_values 和 p_values
            f_values{end+1} = f;
            p_values{end+1} = p;

            % 绘制每个卫星的原始数据
            if iqr
                % 绘制灰色版本以显示 IQR
                plot(f, p, 'Color', [0.7, 0.7, 0.7], 'DisplayName', ['Satellite ' num2str(sat) ' - ' system_band_combined]);
                hold on;

                % 绘制20%阈值线
                noise_threshold = max(p) * 0.2;
                yline(noise_threshold, '--r', 'HandleVisibility', 'off');  % 隐藏后续的图例标签
             
            else
                % 如果未请求 IQR，绘制彩色的原始数据
                plot(f, p, 'DisplayName', ['Satellite ' num2str(sat) ' - ' system_band_combined]);
                hold on;
            end
        end
    end

    % 如果请求 IQR 计算
    if iqr
        % 调用 calculate_weighted_iqr 函数来计算联合功率谱和 IQR
        [combined_left_value, combined_right_value, f_weighted, p_combined, peak_value, peak_freq, sigma,f_combined_valid,consum,consum_total] = calculate_weighted_iqr(f_values, p_values);

        % 在原始图像中应用整体 25%-75% 范围高亮
        for i = 1:length(f_values)
            f = f_values{i};
            p = p_values{i};

            % 在 25%-75% 范围内标记数据
            range_indices = (f >= combined_left_value & f <= combined_right_value);
            plot(f(range_indices), p(range_indices), 'Color', [0.2, 0.6, 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

           % 绘制 IQR 边界线并标注
        xline(combined_left_value, '--b', ['25% Value: ' num2str(combined_left_value, '%.2f')], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left','DisplayName','25% Value');
        xline(combined_right_value, '--r', ['75% Value: ' num2str(combined_right_value, '%.2f')], 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right','DisplayName','75% Value');
        xline(peak_freq, '--k', ['Peak Frequency: ' num2str(peak_freq, '%.2f')], 'LabelVerticalAlignment', 'top','DisplayName','Peak value');
       
    end

    % 设置图例和标题
    legend('show');
    title(plot_title);
    xlabel('Reflector Height [m]');
    ylabel('Volts/Volts');
    grid on;
    hold off;



             
case 'all_data'
    % 所有卫星所有波段数据绘图 (All satellite, all band data plot)
    plot_title = 'All Satellites and Bands Combined';
    f_values = {};   % 存储每个数据集的频率数据
    p_values = {};   % 存储每个数据集的功率数据

    % 获取所有卫星的数据字段 (Get all satellite fields from lsp_results)
    fields = fieldnames(lsp_results);

    % 遍历每个卫星和每个波段的数据 (Iterate over each satellite and each band's data)
    for i = 1:length(fields)
        sat_field = fields{i};
        bands_in_data = fieldnames(lsp_results.(sat_field));

        for j = 1:length(bands_in_data)
            band = bands_in_data{j};
            f = lsp_results.(sat_field).(band).f;
            p = lsp_results.(sat_field).(band).p;

            % 将当前数据集添加到 f_values 和 p_values
            f_values{end+1} = f;
            p_values{end+1} = p;

            % 绘制每个数据集的原始数据
            if iqr
                % 绘制灰色版本以显示 IQR
                plot(f, p, 'Color', [0.7, 0.7, 0.7], 'DisplayName', ['Satellite ' strrep(sat_field, 'sat_', '') ' - ' band]);
                hold on;

                % 绘制20%阈值线
                noise_threshold = max(p) * 0.2;

                yline(noise_threshold, '--g', 'HandleVisibility', 'off');  % 隐藏后续的图例标签

            else
                % 如果未请求 IQR，绘制彩色的原始数据
                plot(f, p, 'DisplayName', ['Satellite ' strrep(sat_field, 'sat_', '') ' - ' band]);
                hold on;
            end
        end
    end

    % 如果请求 IQR 计算
    if iqr
        % 调用 calculate_weighted_iqr 函数来计算联合功率谱和 IQR
        [combined_left_value, combined_right_value, f_weighted, p_combined, peak_value, peak_freq, sigma,f_combined_valid,consum,consum_total] = calculate_weighted_iqr(f_values, p_values);

        % 在原始图像中应用整体 25%-75% 范围高亮
        for i = 1:length(f_values)
            f = f_values{i};
            p = p_values{i};

            % 在 25%-75% 范围内标记数据
            range_indices = (f >= combined_left_value & f <= combined_right_value);
            plot(f(range_indices), p(range_indices), 'Color', [0.2, 0.6, 0.8], 'LineWidth', 1.5,'HandleVisibility', 'off');
        end

        % 绘制 IQR 边界线并标注
        xline(combined_left_value, '--b', ['25% Value: ' num2str(combined_left_value, '%.2f')], 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left','DisplayName','25% Value');
        xline(combined_right_value, '--r', ['75% Value: ' num2str(combined_right_value, '%.2f')], 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right','DisplayName','75% Value');
        xline(peak_freq, '--k', ['Peak Frequency: ' num2str(peak_freq, '%.2f')], 'LabelVerticalAlignment', 'top','DisplayName','Peak value');


    end

    % 设置图例和标题
    legend('show');
    plot(nan, nan, '--g', 'DisplayName', '20% Threshold');  % 创建虚拟绘图以显示图例
    legend('show', 'Location', 'southeast');  % 右下角显示图例
    plot(nan, nan, 'Color', [0.2, 0.6, 0.8], 'LineWidth', 1.5,'DisplayName', '25%-75% Range');
    title(plot_title);
    xlabel('Reflector Height [m]');
    ylabel('Volts/Volts');
    grid on;
    hold off;

    hold off;


    % 最大化图窗 (Maximize the figure window)
    set(figure_handle, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % 确保图像完全展开后再保存 (Ensure the figure is fully expanded before saving)
    drawnow; % 强制完成所有的绘图命令 (Force MATLAB to finish all drawing commands)
            % 新建一个图像窗口并绘制联合功率谱
        
    % Save the figure if the save folder is specified
    if ~isempty(save_folder)
        saveas(gcf, fullfile(save_folder, 'lsp_plot.png'));
    end
    if error_bar
        plot_error_bar(peak_sigma_data, save_folder);
    end

end
    switch plot_type
        case 'single_satellite_band'
           return
        otherwise
        figure;
        plot(f_weighted, p_combined, 'Color', [0, 0.5, 1], 'LineWidth', 1.5, 'DisplayName', 'Combined Power Spectrum');
        hold on;

        % 在联合功率谱图中标注峰值位置和 25%-75% 边界
        % 在联合功率谱图中标注峰值位置和 25%-75% 边界
        xline(peak_freq, '--k', ['Peak Frequency: ' num2str(peak_freq, '%.2f')], ...
               'LabelVerticalAlignment', 'top', 'DisplayName', 'Peak Frequency');
        
        xline(combined_left_value, '--b', ['25% Value: ' num2str(combined_left_value, '%.2f')], ...
               'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'DisplayName', '25% Boundary');
        
        xline(combined_right_value, '--r', ['75% Value: ' num2str(combined_right_value, '%.2f')], ...
               'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'DisplayName', '75% Boundary');

        legend('show');
        title('Combined Power Spectrum');
        xlabel('Reflector Height [m]');
        ylabel('Power');
        grid on;
        % 画出iqr图和左右值 plot iqr and 25%/75% value
        figure 
        y = consum/consum_total;
        plot(f_combined_valid,y)
        hold on
        index1 = find(f_combined_valid == combined_left_value);
        index2 = find(f_combined_valid == combined_right_value);
        yline(y(index1), '--b', ['25% Value: ' num2str(combined_left_value, '%.2f')], ...
               'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'DisplayName', '25% Boundary');
        
        yline(y(index2), '--r', ['75% Value: ' num2str(combined_right_value, '%.2f')], ...
               'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'DisplayName', '75% Boundary');
        title('IQR');
        xlabel('Reflector Height [m]');
        ylabel('Normalized Power');
        grid on;
        hold off
    end
end



function satellite_num = extract_satellite_number(filename)
    % Extracts the satellite number from the filename (assuming a specific format)
    tokens = regexp(filename, 'sat(\d+)', 'tokens');
    satellite_num = str2double(tokens{1});
end

function band = extract_band(filename)
    % Extracts the band information from the filename (assuming a specific format)
    tokens = regexp(filename, '_(L\d)', 'tokens');
    band = tokens{1}{1};
end
function system = extract_system(filename)
    % Extracts the system information (e.g., G, E, R, etc.) from the filename.
    tokens = regexp(filename, '_(L\d)_(\w)', 'tokens');
    system = tokens{1}{2};
end


