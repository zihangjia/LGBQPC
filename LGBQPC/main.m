clear all
warning('off', 'all');
%% 初始化
initialize_path();
load("LGBQPC_parameter_setting.mat");
%% 计算聚类结果
dataset_num = 30;
performance = struct('nmi', zeros(dataset_num, 1), 'ari', zeros(dataset_num, 1), 'time', zeros(dataset_num, 1));
for no = 1:30
    fprintf('计算到第 %d 个数据集\n', no);
    file_name = ['D', num2str(no), '.mat'];
    load(file_name);
    data = (data - min(data, [], 1)) ./ (max(data, [], 1) - min(data, [], 1) + eps);
    tic
    label_pred = GBClustering.lgbqpc(data, class_num, parameter_setting(no, 1), parameter_setting(no, 2));
    performance.time(no) = toc;
    performance.nmi(no) = py.sklearn.metrics.normalized_mutual_info_score(label, label_pred) * 100;
    performance.ari(no) = py.sklearn.metrics.adjusted_rand_score(label, label_pred) * 100;
end
fprintf('%.2f & %.2f & %.3f', mean(max(performance.nmi, [], 2)), mean(max(performance.ari, [], 2)), sum(max(performance.time, [], 2)));
save('experimental_results', 'performance');
%% 函数区
function [] = initialize_path()
% 初始化路径
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, '..', '\datasets'));
end
