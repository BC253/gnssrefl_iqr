function plot_error_bar(peak_sigma_data, save_folder)
    % 获取所有的 band_system_combined 字段
    field_names = fieldnames(peak_sigma_data);

    % 初始化两个图
    % 图1：每个 band_system_combined 的 Error Bar
    figure1 = figure;
    hold on;

    % 图2：系统平均的 Error Bar
    figure2 = figure;
    hold on;

    % 初始化存储系统平均数据的结构
    system_avg_data = struct();

    % 遍历每个 band_system_combined
    for i = 1:length(field_names)
        % 获取当前 band_system_combined 的名称和数据
        band_system_combined = field_names{i};
        data = peak_sigma_data.(band_system_combined);

        % 检查数据是否存在
        if ~isempty(data.satellite_ids)
            % 计算上下误差范围
            error_lower = data.sigma_vals / 2;
            error_upper = data.sigma_vals / 2;

            % 绘制到图1
            figure(figure1);
            errorbar(data.satellite_ids, data.peak_vals, error_lower, error_upper, 'o', ...
                'DisplayName', band_system_combined);

            % 提取系统和波段信息
            split_name = split(band_system_combined, '_');
            band = split_name{1}; % 波段
            system = split_name{2}; % 系统

            % 统计系统平均数据
            key = [band, '_', system];
            if ~isfield(system_avg_data, key)
                system_avg_data.(key) = struct('satellite_ids', [], 'peak_vals', [], 'sigma_vals', []);
            end
            system_avg_data.(key).satellite_ids = [system_avg_data.(key).satellite_ids, data.satellite_ids];
            system_avg_data.(key).peak_vals = [system_avg_data.(key).peak_vals, data.peak_vals];
            system_avg_data.(key).sigma_vals = [system_avg_data.(key).sigma_vals, data.sigma_vals];
        end
    end

    % 绘制系统平均的 Error Bar 到图2
    keys = fieldnames(system_avg_data);
    for j = 1:length(keys)
        avg_key = keys{j};
        avg_data = system_avg_data.(avg_key);

        % 计算平均值和误差
        avg_peak = mean(avg_data.peak_vals);
        avg_sigma = mean(avg_data.sigma_vals);

        % 使用系统 ID 的中点作为 X 轴
        avg_x = mean(avg_data.satellite_ids);

        % 绘制到图2
        figure(figure2);
        errorbar(avg_x, avg_peak, avg_sigma / 2, avg_sigma / 2, 's', 'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [avg_key, '']);
    end

    % 设置图1的图例和标题
    figure(figure1);
    legend('show');
    title('Error Bar for Satellite Groups and Bands');
    xlabel('Satellite ID');
    ylabel('Reflector Height [m]');
    grid on;

    % 设置图2的图例和标题
    figure(figure2);
    legend('show');
    title('Average Error Bar for Satellite Systems and Bands');
    xlabel('Satellite ID');
    ylabel('Reflector Height [m]');
    grid on;

    % 保存图像（如果指定了保存路径）
    if ~isempty(save_folder)
        saveas(figure1, fullfile(save_folder, 'error_bar_individual.png'));
        saveas(figure2, fullfile(save_folder, 'error_bar_avg.png'));
    end
end
