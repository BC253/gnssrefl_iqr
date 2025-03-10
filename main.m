


clear all
clc 
close all


%在这里选择文件夹路径  Select the folder path here
folder_path = 'D:\Amaster\HiWi\GNSSR\bc\refl_code\2025\arcs\ess6\014';
files = dir(fullfile(folder_path, '*.txt'));
filepaths = fullfile({files.folder}, {files.name});
%在这里选择感兴趣的卫星  Select the satellite of interest here 
% 如果后面选单个卫星单个波段，卫星号请一个一个给 If you choose 'single_satellite_band',plz give
% the number one by one
satellite_numbers = [1:1:400];
%在这里选择波段 Select the band of interest here
% bands = {'L1'};
bands = {'L1','L2','L5','L6','L7'};
%在这里选择想要的plot方式，目前有四种
%Select the desired plotting method here, currently there are four Plot :S
%'single_satellite_band', 
%'same_satellite_different_bands',
%'different_satellites_same_band', 
%'all_data'
plot_type = 'all_data';

% 指定保存图片的文件夹路径，若不需要保存则将 save_folder 设为 ''
% Specify the path of the folder where you want to save the image, if you don't need to save it, 
% set save_folder to ''.

save_folder = '';

%use or no use iqr
iqr = true;

error_bar = true; % 设置是否绘制 use or no use Error Bar

plot_gnss_data(filepaths, satellite_numbers, bands, plot_type, save_folder, iqr, error_bar);

%% results

station_ids = [ 2, 3, 4, 5,6,7,8]; 
peak_values = [3.15, 3.01, 2.57, 3.52, 3.35,2.62,2.61]; 
sigma_values = [0.07, 0.11, 0.32, 0.03, 0.04,0.76,0.04];


figure;
hold on;
for i = 1:length(station_ids)
    errorbar(station_ids(i), peak_values(i), sigma_values(i) / 2, sigma_values(i) / 2, 'o', ...
        'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', ['Station ', num2str(station_ids(i))]);
end

xticks(1:length(station_ids)); 
xticklabels(compose('Station %d', station_ids)); 
legend('show');
title('Error Bar for Stations in Motor-Yacht-Club 14.01.25 ');
xlabel('Stations');
ylabel('Peak Height');
grid on;

station_ids = [ 2, 3, 4, 5,6,7,8]; 
peak_values = [0.9, 1.99, 2.58, 3.47, 3.35,2.8,1.76]; 
sigma_values = [0.94, 1.09, 0.08, 1.15, 1.25,0.05,0.58];


figure;
hold on;
for i = 1:length(station_ids)
    errorbar(station_ids(i), peak_values(i), sigma_values(i) / 2, sigma_values(i) / 2, 'o', ...
        'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', ['Station ', num2str(station_ids(i))]);
end

xticks(1:length(station_ids)); 
xticklabels(compose('Station %d', station_ids)); 
legend('show');
title('Error Bar for Stations near bridge 14.01.25 ');
xlabel('Stations');
ylabel('Peak Height');
grid on;

