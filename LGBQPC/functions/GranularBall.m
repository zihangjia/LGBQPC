classdef GranularBall
    %GRANULARBALL 用于实例化粒球并集成与粒球相关的操作
    %
    %  method 1: (Constructor)
    %  obj = GranularBall(in_dataset, in_instance_index)
    %  输入数据集 (in_dataset) 和样本在原始数据集的索引 (in_instance_index),
    %  可以实例化 1 个粒球
    %
    %  method 2:
    %  [gb1, gb2] = self_divide(obj)
    %  将现有粒球 (obj) 划分为 2 个新的子粒球 (gb1, gb2)
    %
    %  method 3:
    %  [gb1, gb2] = self_divide_2_means(obj)
    %  利用 2-均值聚类将现有粒球 (obj) 划分为 2 个新的子粒球 (gb1, gb2)
    %
    %  method 4:
    %  [obj, gb1, gb2] = establish_relationship(obj, gb1, gb2)
    %  父粒球 (obj) 调用本方法, 输入它的 2 个子粒球 (gb1, gb2),
    %  更新 3 者之间的关系并将更新关系之后的 obj, gb1, gb2 返回
    %
    %  method 5: out_gamma = update_granularity_level(obj, gb1, gb2, in_lambda, in_gamma)
    %  父粒球 (obj) 调用本方法, 输入它的 2 个子粒球 (gb1, gb2),
    %  惩罚系数 (in_lambda), 当前粒度水平范围 (in_gamma)
    %  根据构造父粒球的数据集的特征的相关系数和特征方差的变异系数决定是否一定要分裂这个粒球
    %  如果一定要分裂这个粒球, 那么更新粒度水平的取值范围

    properties
        % 构成粒球的基本属性
        dataset            % 构造粒球的样本构成的矩阵 (n*m, n 个样本, m 个特征)
        instance_index     % 构造粒球的样本在原始数据集中的索引 (n*1, n 个样本对应 n 个索引)
        instance_num       % 构造粒球样本的个数 (1*1)
        feature_num        % 构造粒球特征的个数

        % 粒球的基本属性
        center             % 粒球的球心 (1*m, 球心可以看作具有 m 个特征的样本)
        radius_max         % 粒球的最大半径 (1*1)
        radius_ave         % 粒球的平均半径 (1*1)

        % 通过计算得到的高级属性
        coverage           % 与球心的距离小于平均半径的样本的数量 (1*1)
        granularity_level  % 粒球的粒度水平 (1*1)
        specificity        % 粒球的特异性 (1*1, $f_{2}(R_{X}^{ave})$)
        quality            % 粒球的质量 (1*1, coverage * specificity - lambda）

        % 计算过程中产生的属性
        distance_vector    % 每个样本与球心的距离 (n*1)

        % 粒球与其它粒球之间的关系
        is_leaf            % 当前粒球是否是叶节点 (1*1, boolean)
        num                % 当前粒球编号 (1*1)
        num_father         % 父粒球的编号 (1*1)
        num_left_child     % 左孩子节点编号 (1*1)
        num_right_child    % 右孩子节点编号 (1*1)


        % **********最佳质量和最佳子粒球组合*********
        best_quality       % 粒球的最佳质量 (1*1)
        best_combination   % 粒球的最佳子粒球组合 (t*1, t 是最佳子粒球个数)

        % ******************** 标签 ********************
        label_pred         % 预测的标签
    end

    methods
        function obj = GranularBall(in_dataset, in_instance_index)
            %GRANULARBALL 构造此类的实例
            % param: in_dataset 构造粒球的数据集 (n*m)
            % param: in_instance_index 构造粒球的样本在原始数据集的索引 (n*1)
            % out: obj 生成的粒球

            % ********************** 默认初始化 **********************
            if nargin == 0
                return
            end
            % 如果没有指定输入的数据集, 此时是为了初始化占据空间
            % ********************** 默认初始化 **********************

            % *********** 初始化构造粒球的基本信息 ***********
            obj.dataset = in_dataset;
            [obj.instance_num, obj.feature_num] = size(obj.dataset);
            % (1) 初始化构造粒球的数据集 (obj.dataset)
            % (2) 获取数据集的样本个数 (obj.instance_num) 和特征个数 (obj.feature_num)

            if nargin == 1
                obj.instance_index = (1:obj.instance_num)';
            else
                obj.instance_index = in_instance_index;
            end
            % (1) 如果缺省 input_instance_index,
            %     那么默认是输入原始数据集构造粒球, 因此用 (1:obj.instance_num)' 作为样本索引
            % (2) 如果没有缺省 input_instance_index,
            %     那么用输入的样本索引初始化粒球的样本索引
            % *********** 初始化构造粒球的基本信息 ***********

            % *********** 计算粒球的基本属性 ***********
            obj.center = mean(obj.dataset, 1);
            obj.distance_vector = sqrt(sum((obj.dataset - obj.center) .^ 2, 2));
            obj.radius_max = max(obj.distance_vector);
            obj.radius_ave = mean(obj.distance_vector);
            % (1) 球心是样本在各个特征下的均值构成的向量
            % (2) 距离向量是样本到球心的距离
            % (3) 最大半径是样本到球心距离的最大值
            % (4) 平均半径是样本到球心距离的均值
            % *********** 计算粒球的基本属性 ***********

            % *********** 初始化粒球与其它粒球的关系 ***********
            obj.num = 1;
            obj.num_father = 0;
            obj.num_left_child = -1;
            obj.num_right_child = -1;
            obj.is_leaf = true;
            % (1) 默认设置: 粒球编号为 1
            % (2) 默认设置: 粒球父节点编号为 0
            % (3) 默认设置: 粒球的左孩子节点为 -1
            % (4) 默认设置: 粒球的右孩子节点为 -1
            % (5) 默认设置: 粒球在粒球二叉树中是叶节点
            % *********** 初始化粒球与其它粒球的关系 ***********

            % *********** 初始化粒球的高级属性 ***********
            obj.coverage = sum(obj.distance_vector <= obj.radius_ave);
            % (1) 粒球的覆盖范围是和球心距离小于平均半径的样本数量

            % *********** 初始化预测的标签 **************
            obj.label_pred = -1;
        end

        function [gb1, gb2] = self_divide(obj)
            % 当前粒球 (obj) 划分为 2 个粒球 (gb1, gb2)
            % param: obj 当前粒球
            % out: [gb1, gb2] 2 个新粒球构成的向量

            % *********** 计算 point_1 和 point_2 ***********
            % point_1 是离球心最远的样本, point_2 是离 point_1 最远的样本

            [~, point_1_index] = max(obj.distance_vector);
            point_1 = obj.dataset(point_1_index, :);
            distance_point_1 = sqrt(sum((obj.dataset - point_1) .^ 2, 2));
            % (1) 计算 point_1 的索引 (point_1_index), point_1 是离球心最远的样本
            % (2) 得到 point_1 对应的行向量 (point_1)
            % (3) 计算每个样本与 point_1 的距离 (dinstance_point_1)

            [~, point_2_index] = max(distance_point_1);
            point_2 = obj.dataset(point_2_index, :);
            distance_point_2 = sqrt(sum((obj.dataset - point_2) .^ 2, 2));
            % (1) 计算 point_2 的索引 (point_2_index), point_2 是离 point_1 最远的样本
            % (2) 得到 point_2 对应的行向量 (point_2)
            % (3) 计算每个样本与 point_2 的距离 (dinstance_point_2)

            nearest_point_1 = distance_point_1 < distance_point_2;
            nearest_point_2 = ~nearest_point_1;
            % (1) 获取离 point_1 更近的样本的逻辑向量 (nearest_point_1)
            % (2) 获取离 point_2 更近的样本的逻辑向量 (nearest_point_1)
            gb1 = GranularBall(obj.dataset(nearest_point_1, :), obj.instance_index(nearest_point_1));
            gb2 = GranularBall(obj.dataset(nearest_point_2, :), obj.instance_index(nearest_point_2));

            % (1) 实例化子球 gb1, 由离 point_1 更近的样本构成
            % (2) 实例化子球 gb2, 由离 point_2 更近的样本构成
        end

        function [gb1, gb2] = self_divide_2_means(obj, initialization)
            % 当前粒球 (obj) 划分为 2 个粒球 (gb1, gb2)
            % param: obj 当前粒球
            % out: [gb1, gb2] 2 个新粒球构成的向量

            if nargin < 2 % 如果没有指定 2-均值聚类的初始化方法, 则采用固定的方法初始化
                initialization = 'fixed';
            end

            if strcmp(initialization, 'random')
                label = kmeans(obj.dataset, 2, 'Replicates', 3);
                % 随机挑选 2-均值聚类的初始点, 重复 3 次取目标函数值最低的
            else
                % *********** 计算 point_1 和 point_2 ***********
                % point_1 是离球心最远的样本, point_2 是离 point_1 最远的样本

                [~, point_1_index] = max(obj.distance_vector);
                point_1 = obj.dataset(point_1_index, :);
                distance_point_1 = sqrt(sum((obj.dataset - point_1) .^ 2, 2));
                % (1) 计算 point_1 的索引 (point_1_index), point_1 是离球心最远的样本
                % (2) 得到 point_1 对应的行向量 (point_1)
                % (3) 计算每个样本与 point_1 的距离 (dinstance_point_1)

                [~, point_2_index] = max(distance_point_1);
                point_2 = obj.dataset(point_2_index, :);
                % (1) 计算 point_2 的索引 (point_2_index), point_2 是离 point_1 最远的样本
                % (2) 得到 point_2 对应的行向量 (point_2)
                % (3) 计算每个样本与 point_2 的距离 (dinstance_point_2)

                mid_point = (obj.center + [point_1; point_2]) / 2;
                label = kmeans(obj.dataset, 2, 'Start', mid_point);
                % 以球心和 point_1 与 point_2 的中点作为 2-均值聚类的七十点
            end
            gb1 = GranularBall(obj.dataset(label == 1, :), obj.instance_index(label == 1));
            gb2 = GranularBall(obj.dataset(label == 2, :), obj.instance_index(label == 2));
            % (1) 实例化子球 gb1, 由离 point_1 更近的样本构成
            % (2) 实例化子球 gb2, 由离 point_2 更近的样本构成
        end

        function [obj, gb1, gb2] = establish_relationship(obj, gb1, gb2)
            % 更新当前粒球的与子粒球的关系
            % param: input_gb1 子粒球1, 设置为左孩子节点
            % param: input_gb2 子粒球2, 设置为右孩子节点
            % out: [obj, gb1, gb2] 父粒球, 子粒球1, 子粒球2 构成的向量

            % ********* 更新粒球在二叉树中的角色 *********
            obj.is_leaf = false;
            gb1.is_leaf = true;
            gb2.is_leaf = true;
            % (1) 将父粒球 (obj) 设置为非叶节点
            % (2) 将子粒球1 (gb1) 设置为叶节点
            % (3) 将子粒球2 (gb2) 设置为叶节点

            % ********* 更新粒球之间的关系 *********
            obj.num_left_child = gb1.num;
            obj.num_right_child = gb2.num;
            gb1.num_father = obj.num;
            gb2.num_father = obj.num;
            % (1) 将父粒球 (obj) 左孩子节点编号设置为子粒球1的编号 (gb1.num)
            % (2) 将父粒球 (obj) 右孩子节点编号设置为子粒球2的编号 (gb2.num)
            % (3) 将子粒球1的父节点编号设置为父粒球编号 (obj.num)
            % (4) 将子粒球2的父节点编号设置为父粒球编号 (obj.num)
        end

        function out_gamma = update_granularity_level(obj, gb1, gb2, in_penalty_coefficient, in_gamma, significance_level)
            % 根据当前粒球 (obj) 和其子粒球 (gb1, gb2) 来更新粒度水平的取值范围 (out_gamma)
            % param: gb1 子粒球 1
            % param: gb2 子粒球 2
            % param: in_lambda 粒球质量的惩罚系数
            % param: in_gamma 原来的粒度水平取值范围
            % out: out_gamma 更新之后的粒度水平取值范围

            is_not_spherical = Tool.bartlett_sphericity_test(obj.dataset, significance_level);

            % ********************* 计算粒度水平的值 *********************
            if is_not_spherical
                % 如果当前粒球的数据集不是球形的
                % 此时可能需要更新 gamma 的取值范围

                % ************* 计算 3 次不等式的系数 *************
                % 令 f(x) = ax^3 + bx^2 + cx + d
                a = -1 * in_penalty_coefficient * obj.radius_ave * gb1.radius_ave * gb2.radius_ave;
                b = obj.radius_ave * (gb1.coverage+gb2.coverage) + gb1.radius_ave * (gb2.coverage-obj.coverage) + gb2.radius_ave * (gb1.coverage-obj.coverage) - in_penalty_coefficient * (obj.radius_ave+gb1.radius_ave+gb2.radius_ave);
                c = obj.radius_ave * (gb1.coverage+gb2.coverage-in_penalty_coefficient) + gb1.radius_ave * (gb2.coverage-obj.coverage-in_penalty_coefficient) + gb2.radius_ave * (gb1.coverage-obj.coverage-in_penalty_coefficient);
                d = gb1.coverage + gb2.coverage - obj.coverage - in_penalty_coefficient;
                % (1)-(4) 计算 3 次不等式的各项系数

                solution = Tool.solve_3_order_polynomial_inequality([a b c d]);
                % (1) 通过 Tool.solve_3_order_polynomial_inequality([a b c d]) 函数求 3
                %     次不等式的解空间, 0, 1, 2 次不等式是解析解, 3 次不等式是数值解, f(x) > 0

                out_gamma = Tool.interval_intersection(solution, in_gamma);
                % (1) 根据 3 次不等式的解空间更新粒度水平的取值范围, 新的取值范
                %     围是解空间和原取值范围的交集

                if size(out_gamma, 1) == 0
                    out_gamma = in_gamma;
                end
                % (1) 如果更新之后的粒度水平取值范围为空
                %     那么不更新粒度水平了, 即输出的粒度水平范围为输入的粒度水平范围
            else
                out_gamma = in_gamma;
                % (1) 此时无需更新 gamma, 即 out_gamma = input_gamma
            end
        end

        
    end
end

