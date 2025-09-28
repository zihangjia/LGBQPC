classdef GBClustering
    methods(Static)
        function [label_pred, gbs] = lgbqpc(in_feature, in_class_num, in_penalty_coefficient, in_neighbor_num, in_significance_level)
            % 实现 Local Granular-Ball Quality Peaks Clutering Based on
            % k-nearest neighbor graph
            % param: in_feature:                   样本的特征矩阵
            % param: in_class_num:                 样本的真实类别数
            % param: in_penalty_coefficient:       粒球生成模型的惩罚系数
            % param: in_neighbor_num:              k-近邻图邻居的数量
            % param: in_significance_level:        球形检验的显著性水平

            % ********************************* 初始化 *********************************
            [instance_num, ~] = size(in_feature);
            label_pred = zeros(instance_num, 1);
            if nargin < 3 || isempty(in_penalty_coefficient)
                in_penalty_coefficient = 0.01 * power(instance_num, 1/3);
            end
            if nargin < 4 || isempty(in_neighbor_num)
                in_neighbor_num = [];
            end
            if nargin < 5 || isempty(in_significance_level)
                in_significance_level = 0.01;
            end
            % (1) 用 instance_num 存储样本个数
            % (2) 用 label_pred 存储每个样本的标签, 初始化为 instance_num * 1
            %     的向量
            % (3)-(5) 如果第 3 个参数缺省, 则设置粒球生成模型的惩罚系数为
            %         0.01*(样本数量开 3 次方)
            % (6)-(8) 如果第 4 个参数缺省, 则 k-近邻图邻居的数量设置为 1
            % ********************************* 初始化 *********************************


            gbs = GBClustering.generate_granular_balls(in_feature, in_penalty_coefficient, in_significance_level);
            clear in_feature
            % 根据样本特征和惩罚系数生成粒球并存储粒球相关的信息

            gbs_neighbor_graph = GBClustering.generate_knn_graph(gbs.gbs_distance, in_neighbor_num);
            % 根据粒球之间的距离和粒球的 k 近邻建立粒球邻域图

            gbs = GBClustering.calculate_gbs_relative_quality(gbs, gbs_neighbor_graph);
            % 计算粒球的相对质量

            gbs = GBClustering.calculate_gbs_relative_distance(gbs, gbs_neighbor_graph);
            % 计算粒球的相对距离

            label_pred = GBClustering.obtain_instance_label(gbs, in_class_num, label_pred);
            % 计算样本的标签
        end

        function gbs = generate_granular_balls(in_feature, in_lambda, in_significance_level)
            % 根据样本的特征矩阵 (in_feature) 和惩罚系数 (in_lambda) 生成粒球,
            % 返回与粒球相关的信息构成的结构体 (gbs)
            % param: in_feature: 样本的特征矩阵
            % param: in_lambda:  生成粒球的惩罚系数
            % out: gbs 存储粒球信息的结构体

            gbs.content = GBComputing.GB_POJG_PLUS(in_feature, in_lambda, in_significance_level);
            gbs.num = size(gbs.content, 1);
            gbs.center = reshape([gbs.content.center], [], gbs.num)';
            gbs.radius_ave = [gbs.content.radius_ave]';
            gbs.center_distance = pdist2(gbs.center, gbs.center, 'fasteuclidean') ;
            gbs.gbs_distance = max(1e-8 * ones(gbs.num), gbs.center_distance - (gbs.radius_ave + gbs.radius_ave'));
            gbs.quality = [gbs.content.quality]' + in_lambda;
            gbs.instance_index = {gbs.content.instance_index}';
            clear gbs.content
            % (1) 根据样本的特征矩阵 (in_feature) 和惩罚系数 (in_lambda) 生成
            %     粒球, 存储于 gbs.content
            % (2) 用 gbs.num 存储生成的粒球的个数
            % (3) 用 gbs.center 存储生成的粒球的球心 (gbs.num * feature_num):
            %     用 [gbs.center] 获得由所有粒球球心构成的行向量, 然后利用
            %     reshape 函数, 生成 gbs.num 列的矩阵, 每一列都对应一个粒球的
            %     球心, 再通过转置, 获得矩阵 gbs.center. gbs.center(i, j) 表
            %     示第 i 个粒球的球心的第 j 个分量
            % (4) 用 gbs.radius_ave 存储生成的粒球的平均半径 (gbs.num * 1)
            % (5)-(6) 用 gbs.gbs_distance 存储粒球之间的距离 (gbs.num * gbs.num)
            % (7) 用 gbs.quality 存储粒球的质量 (gbs.num * 1), 把粒球生成模型的惩罚系数加回去
            % (8) 存储构造每个粒球对应的样本标签索引构成的元胞数组
            % (9) 释放生成的粒球及对应的数据集
        end

        function knn_graph = generate_knn_graph(in_distance_matrix, in_neighbor_num)
            % 根据距离矩阵和邻居数量建立 k-近邻图
            % param: in_distance_matrix 距离矩阵
            % param: in_neighbor_num 邻居数量
            % out: knn_graph k-近邻图

            node_num = size(in_distance_matrix, 1);
            if isempty(in_neighbor_num)
                in_neighbor_num = floor(log2(node_num)) + 1;
            end
            in_distance_matrix(1:node_num+1:end) = inf;
            if  in_neighbor_num >= node_num - 1
                knn_graph = graph(in_distance_matrix, 'upper');
                return
            end
            % (1)     获取节点数量
            % (2)-(4) 如果未指定最近邻个数, 则设置为 log2(gbs.num)+1 向下取整
            % (5)     将自身距离设置为无穷大，防止被选为邻居, 对角线索引的差为 node_num + 1
            % (6)-(9) 如果邻居数量大于等于节点数量-1，则直接返回全连接图,
            %         使用距离矩阵的上三角阵生成图，避免重复边

            [nearest_distance, nearest_neighbor] = mink(in_distance_matrix, in_neighbor_num, 2);
            edge.begin  = repmat((1:node_num)', in_neighbor_num, 1);  % 边的起点列向量
            edge.end    = nearest_neighbor(:);                        % 边的终点列向量
            edge.weight = nearest_distance(:);                        % 边的权重列向量
            adj_matrix = sparse(edge.begin, edge.end, edge.weight, node_num, node_num);
            adj_matrix = max(adj_matrix, adj_matrix');
            knn_graph = graph(adj_matrix, 'upper');
            % (1) 使用 mink 函数取每个节点最近的 in_neighbor_num 个邻居,
            %     nearest_distance 存储每个节点最近的 in_neighbor_num 个邻居的距离
            %     nearest_neighbor 存储每个节点最近的 in_neighbor_num 个邻居的编号
            % (2)-(4) 向量化构造边
            % (5) 构造稀疏邻接矩阵
            % (6) 确保邻接矩阵对称
            % (7) 根据邻接矩阵构造 k-近邻图
        end


        function gbs = calculate_gbs_relative_quality(gbs, gbs_neighbor_graph)
            % 输入粒球信息和粒球邻域图, 输出粒球的相对质量向量
            % param: gbs: 结构体, 存储关于粒球的各种信息
            % param: gbs_neighbor_graph: 粒球邻域图

            gbs_adjacency_matrix = logical(adjacency(gbs_neighbor_graph));
            neighbor_sum_quality = gbs_adjacency_matrix * gbs.quality(:);
            neighbor_count = sum(gbs_adjacency_matrix, 2);
            gbs_neighbor_mean_quality = neighbor_sum_quality ./ (neighbor_count + eps);
            gbs.relative_quality = gbs.quality(:) ./ (gbs_neighbor_mean_quality + eps);
            normalized_quality = (gbs.quality - min(gbs.quality)) ./ (max(gbs.quality) - min(gbs.quality) + eps);
            normalized_relative_quality = (gbs.relative_quality - min(gbs.relative_quality)) ./ (max(gbs.relative_quality) - min(gbs.relative_quality) + eps);
            gbs.combined_quality = normalized_quality + normalized_relative_quality;
            % (1) 获取稀疏邻接矩阵, 数据类型为逻辑型
            % (2) 计算每个粒球的邻居的质量和
            % (3) 计算每个粒球的邻居的数量
            % (4) 计算每个粒球的邻居的平均质量, eps 防止除 0
            % (5) 计算每个粒球的相对质量
            % (6) 归一化粒球的质量
            % (7) 归一化粒球的相对质量
            % (8) 计算每个粒球的组合质量
            % (9) 将粒球的组合质量归一化到 [0,1] 区间
        end


        function gbs = calculate_gbs_relative_distance(gbs, gbs_neighbor_graph)
            % 根据粒球邻域图的测地距离和粒球组合质量计算粒球之间的相对距离

            gbs.geodesic_distance = distances(gbs_neighbor_graph);
            [~, combined_quality_index] = sort(gbs.combined_quality, 'descend');
            gbs_relative_distance = zeros(gbs.num, 1);
            gbs_nearest_neighbor = -1 * ones(gbs.num, 1);
            % (1) 由粒球邻域图得到粒球之间的测地距离
            % (2) 对粒球的组合质量降序排列, combined_quality_index(i) 存储第 i 大的粒球索引
            % (3) 初始化粒球相对距离向量
            % (4) 初始化粒球的最近邻

            for item = 2:gbs.num
                gb_now_index = combined_quality_index(item);
                gbs_higher_quality_index = combined_quality_index(1:item-1);
                [gbs_relative_distance(gb_now_index), idx] = min(gbs.geodesic_distance(gb_now_index, gbs_higher_quality_index));
                gbs_nearest_neighbor(gb_now_index) = gbs_higher_quality_index(idx);
                % (1) 当前粒球索引
                % (2) 比当前粒球组合质量高的粒球索引
                % (3) 得到当前粒球的相对距离和相对最近邻在比他组合质量大的粒球的索引
                % (4) 存储当前粒球的最近邻
            end
            gbs_relative_distance(combined_quality_index(1)) = max(gbs_relative_distance);
            gbs_nearest_neighbor(combined_quality_index(1)) = combined_quality_index(1);
            % (1)-(2) 组合质量最大的粒球的相对距离是最大的粒球相对距离, 最近邻是自身

            % 归一化相对距离
            gbs.relative_distance = gbs_relative_distance;
            gbs.combined_quality_index = combined_quality_index;
            gbs.nearest_neighbor = gbs_nearest_neighbor;
            % (3) gbs.relative_distance: 粒球的相对距离
            % (4) gbs.nearest_neighbor: 粒球的最近邻 (组合质量比它高且距离最近)
            % (5) gbs.combined_quality_index: 粒球按组合质量降序排列的索引
        end


        function label_pred = obtain_instance_label(gbs, in_class_num, label_pred)
            % 基于粒球相对质量和相对距离为粒球及其样本分配标签
            % 完全保留原有逻辑，最后一部分向量化
            % param gbs: 粒球结构体
            % param in_class_num: 类别数
            % param label_pred: 样本标签向量 (预分配)
            % out label_pred: 更新后的样本标签

            decision_value = gbs.combined_quality .* gbs.relative_distance;
            [~, decision_index] = sort(decision_value, 'descend');
            gbs_label = zeros(gbs.num, 1);
            % (1) 基于粒球的相对质量和相对距离计算决策值
            % (2) 对决策值从高到低进行排序, decision_index(i) 表示决策值排名为 i 的粒球的索引
            % (3) 初始化粒球的质量标签

            for item = 1:min(in_class_num, gbs.num)
                gbs_label(decision_index(item)) = item;
            end
            % 将决策值排名为 1-in_class_num 的粒球赋予标签

            % 按 combined_quality_index 顺序遍历粒球，确保传播标签时最近邻已赋值
            for item = 1:gbs.num
                idx = gbs.combined_quality_index(item);
                if gbs_label(idx) == 0
                    gbs_label(idx) = gbs_label(gbs.nearest_neighbor(idx));
                end
            end
            % 遍历每个粒球，如果该粒球没有标签，则其标签与最近邻一致

            % ===== 向量化赋值样本标签 =====
            instance_counts = cellfun(@numel, gbs.instance_index);
            label_pred(vertcat(gbs.instance_index{:})) = repelem(gbs_label, instance_counts)';
            % (1) 统计每个粒球包含样本数量
            % (2) 将粒球标签重复对应样本数量
            % (3) 直接赋值给样本标签向量 label_pred
        end
    end
end




