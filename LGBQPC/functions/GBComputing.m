classdef GBComputing
    % GBComputing 粒球计算类, 实现与粒球计算有关的静态方法
    % GB_POJG_PLUS 由 5 个子函数构成
    % 1. establish_gbs_tree:             建立粒球二叉树并获得粒度水平的取值范围
    % 2. detect_abnormal_leaf_nodes:     检测粒球二叉树的异常叶节点并分裂异常叶节点
    % 3. update_tree:                    计算粒球二叉树中每个粒球的质量
    % 4. calculate_gbs_best_combination: 计算最佳粒球组合
    % 5. detect_abnormal_gbs:            检测并分裂异常的粒球

    methods(Static)

        function gbs = GB_POJG_PLUS(data_matrix, penalty_coefficient, significance_level)
            % 输入样本的特征矩阵, 惩罚系数和显著性水平生成用于表示样本的粒球
            % param: data_matrix, instance_num * feature_num, 样本的特征矩阵
            % param: penalty_coefficient 粒球计算模型中关于粒球个数的惩罚系数

            % ************************ 参数缺省处理 ************************
            if nargin < 2 || isempty(penalty_coefficient)
                penalty_coefficient = 0.1 * power(size(data_matrix, 1), 1/3);
            end
            if nargin < 3 || isempty(significance_level)
                significance_level = 0.01;
            end
            % 如果第 2 个参数惩罚系数 penalty_coefficient 缺省, 则其默认值为 0.1 * (样本数量开三次方)
            % 如果第 3 个参数显著性水平 significance_level 缺省, 则其默认值维 0.01
            % ************************ 参数缺省处理 ************************

            % ********** 根据样本的特征矩阵和惩罚系数建立粒球二叉树 **********
            [tree, leaf_queue, gamma_range] = GBComputing.establish_gbs_tree...
                (data_matrix, penalty_coefficient, significance_level);
            % (1) 建立粒球二叉树 (tree)
            % (2) 获得存储粒球二叉树叶节点的顺序循环队列 (leaf_queue)
            % (3) 获得粒度水平的取值范围 (gamma_range)
            % ********** 根据样本的特征矩阵和惩罚系数建立粒球二叉树 **********

            % *********** 检测异常叶节点并进一步划分异常叶节点建树 ***********
            [tree, leaf_queue] = GBComputing.detect_abnormal_leaf_nodes(tree, leaf_queue);

            % ************ 根据粒度水平范围更新粒球二叉树 ************
            tree = GBComputing.update_tree(tree, gamma_range, penalty_coefficient);
            % 根据粒度水平范围 (gamma_range) 更新粒球二叉树 (tree) 的质量 (quality)
            % 并初始化最佳质量 (best_quality) 和最佳子粒球组合 (best_combination)

            % ************* 获取最佳粒球组合 *************
            tree = GBComputing.calculate_gbs_best_combination(tree, leaf_queue);
            % (1) 计算粒球二叉树根节点的最佳子粒球组合 (粒球编号构成的数组)
            % (2) 获取根节点的最佳子粒球组合 (粒球数组)

            gbs = GBComputing.detect_abnormal_gbs(tree, penalty_coefficient);
            % 检测异常粒球并分裂异常粒球
        end

        function [tree, leaf_queue, gamma_range] = establish_gbs_tree(data_matrix, penalty_coefficient, significance_level)
            % 输入样本的特征矩阵和惩罚系数返回粒球二叉树, 叶节点队列和粒度水平取值范围
            % param: data_matrix, instance_num * feature_num, 样本的特征矩阵
            % param: penalty_coefficient, 1*1, 惩罚系数
            % out: tree 粒球二叉树 (结构体, 字段 node 存储粒球数组,
            %      字段 node_num 存储节点个数)
            % out: leaf_queue 叶节点队列 (顺序循环队列, 存储粒球二叉树叶节点的编号)
            % out: gamma_range, t*2, 粒度水平的取值范围

            % ********************* 初始化 ********************
            [instance_num, ~] = size(data_matrix);
            tree = struct();
            tree.node = createArray(instance_num * 2, 1, "GranularBall");
            tree.node(1) = GranularBall(data_matrix);
            tree.node_num = 1;
            queue = Queue(instance_num * 2, 'uint32');
            leaf_queue = Queue(instance_num, 'uint32');
            queue = queue.append(1);
            gamma_range = [0, inf];
            % (1) instance_num 样本个数
            % (2) tree 是粒球二叉树, 以结构体的形式实现
            % (3) tree.node 以粒球数组模拟粒球二叉树, createArray(size_1, size_2,
            %     data_type) 创建 size_1 * size_2 大小, 数据类型为 data_type 的矩阵
            % (4) tree.node(1) 粒球二叉树的根节点是原始数据集生成的粒球
            % (5) tree.node_num 存储粒球二叉树的节点个数, 设置为 1
            % (6) queue 是顺序循环队列, 存储仍需划分的粒球节点在粒球二叉树的编号,
            %     这里用数组模拟顺序循环队列, 数组的数据类型为 uint32
            % (7) leaf_queue 是顺序循环队列, 存储粒球二叉树的叶节点在粒球二叉树
            %     的编号, 这里用数组模拟顺序循环队列, 数组的数据类型为 uint32
            % (8) 将粒球二叉树根节点的编号 (1) 加入顺序循环队列 (queue)
            % (9) 初始化粒度水平参数 gamma 的取值范围为 [0, inf]
            % ********************* 初始化 ********************

            % ************ 建立粒球二叉树 ************
            while ~queue.is_empty()
                % 队列 queue 非空时, 说明仍存在粒球过大, 需要被进一步划分, 进而循环执行以下操作

                % ******* 预备工作: 根据队列从粒球树中选定要处理的粒球 *******
                num_now = queue.get();
                queue = queue.delete();
                gb_now = tree.node(num_now);
                % (1) 获取队列 (queue) 队首的粒球编号 (now_num)
                % (2) 删除队列 (queue) 队首的粒球编号
                % (3) 在粒球二叉树 (tree) 中获取节点编号为 now_num 的粒球 (gb_now)
                % ******* 预备工作: 根据队列从粒球树中选定要处理的粒球 *******

                % ********* 更新粒球二叉树或更新粒球二叉树叶节点队列 *********
                if  gb_now.instance_num > power(instance_num, 1/3)
                    % 当构成当前粒球的样本数量大于 (样本数量开三次方) 且有不同的样本, 则当前粒球需要分裂

                    % **************** 更新粒球二叉树和队列 ****************
                    [gb1, gb2] = gb_now.self_divide_2_means();
                    gb1.num = tree.node_num + 1;
                    gb2.num = tree.node_num + 2;
                    [tree.node(num_now), tree.node(gb1.num), tree.node(gb2.num)] = gb_now.establish_relationship(gb1, gb2);
                    tree.node_num = tree.node_num + 2;
                    queue = queue.append(gb1.num);
                    queue = queue.append(gb2.num);
                    % (1) 将当前粒球 (gb_now) 分裂为 2 个子粒球 (gb1, gb2)
                    % (2)-(3) 设置当前粒球 (gb_now) 的 2 个子粒球 (gb1, gb2) 的编号
                    % (4) 建立 gb_now 与 gb1 与 gb2 之间的父子关系并在粒球二叉树 (tree) 中更新对应的节点
                    % (5) 更新粒球二叉树中节点的数量 (tree.node_num)
                    % (6)-(7) 将新生成的粒球的编号 (gb1.num, gb2.num) 纳入队列 queue 中
                    % **************** 更新粒球二叉树和队列 ****************

                    % *************** 更新粒度水平的取值范围 ***************
                    if gb_now.instance_num > power(instance_num, 1/3)
                        gamma_range = gb_now.update_granularity_level(gb1, gb2, penalty_coefficient, gamma_range, significance_level);
                    end
                    % (1) 如果构成当前粒球的样本个数 (gb_now.instance_num) 大于原始数据集样本个数开三次方 (power(instance_num, 1/3))
                    % (2) 那么根据当前粒球 (gb_now), 其生成的子粒球 (gb1, gb2) 和参数 (lambda) 更新粒度水平 (gamma_range)
                    % *************** 更新粒度水平的取值范围 ***************
                else
                    % 当前粒球不需要分裂
                    leaf_queue = leaf_queue.append(num_now);
                    % 将当前节点的编号纳入叶节点队列 leaf_queue
                end
            end
        end

        function [tree, leaf_queue] = detect_abnormal_leaf_nodes(tree, leaf_queue)
            % 检测异常叶节点, 如果一个叶节点是异常叶节点, 那么该叶节点需要被划分
            % 异常粒球的定义 (二者满足一个即可):
            % (1) 粒球的最大半径 > 粒球最大半径的均值 + 粒球最大半径的标准差
            % (2) 粒球的平均半径 > 粒球平均半径的均值 + 粒球平均半径的标准差
            %     且粒球的样本数量 小于 粒球样本数量的均值 - 粒球样本数量的标准差
            % param: tree 粒球二叉树
            % param: leaf_queue 存储叶节点编号的队列

            % ************************* 初始化 *************************
            leaf_node_num = leaf_queue.element(1:leaf_queue.rear-1);
            value.radius_max = [tree.node(leaf_node_num).radius_max]';
            value.radius_ave = [tree.node(leaf_node_num).radius_ave]';
            value.instance_num = [tree.node(leaf_node_num).instance_num]';
            [value.radius_max_std, value.radius_max_mean] = std(value.radius_max);
            [value.radius_ave_std, value.radius_ave_mean] = std(value.radius_ave);
            [value.instance_num_std, value.instance_num_mean] = std(value.instance_num);
            leaf_queue = Queue(tree.node(1).instance_num, "uint32");
            queue = Queue(2 * tree.node(1).instance_num, "uint32");
            % (1) leaf_node_num 存储所有叶节点的编号
            % (2) value.radius_max 存储所有叶节点最大半径构成的列向量
            % (3) value.radius_ave 存储所有叶节点平均半径构成的列向量
            % (4) value.instance_num 存储所有叶节点的样本数量构成的列向量
            % (5)-(7) 计算叶节点最大半径、平均半径以及样本数量的标准差和均值
            % (8) 重新定义叶节点编号的队列, 队列 leaf_queue 最大包含样本个数个样本
            % (9) 初始化异常叶节点编号的队列, 队列 queue 最大包含 2 倍样本个数个样本
            % ************************* 初始化 *************************

            % ************** 在队列中存储异常叶节点和正常叶节点 **************
            abnormal_max_index = value.radius_max > (value.radius_max_mean + value.radius_max_std);
            abnormal_ave_index = (value.radius_ave > (value.radius_ave_mean + value.radius_ave_std)) & ...
                (value.instance_num < (value.instance_num_mean - value.instance_num_std));
            abnormal_index = abnormal_ave_index | abnormal_max_index;
            queue.element = leaf_node_num(abnormal_index);
            queue.now_element_num = sum(abnormal_index);
            queue.rear = queue.now_element_num + 1;
            leaf_queue.element = leaf_node_num(~abnormal_index);
            leaf_queue.now_element_num = sum(~abnormal_index);
            leaf_queue.rear = leaf_queue.now_element_num + 1;
            % (1) - (3) 获取异常叶节点的索引 abnormal_index
            % (4) - (6) 初始化存储异常节点的队列 queue
            % (7) - (9) 初始化存储正常叶节点的队列 leaf_queue
            % ************** 在队列中存储异常叶节点和正常叶节点 **************

            while ~queue.is_empty() % 当异常节点非空时, 循环划分异常节点
                num_now = queue.get();
                gb_now = tree.node(num_now);
                queue = queue.delete();
                % (1) 从队列获取异常节点的编号 num_now
                % (2) 获取当前粒球
                % (3) 从队列中删除队首节点
                if gb_now.radius_max > (value.radius_max_mean + value.radius_max_std) || ...
                        (gb_now.radius_ave > (value.radius_ave_mean + value.radius_ave_std)) && ...
                        (gb_now.instance_num < (value.instance_num_mean - value.instance_num_std))
                    % 如果当前叶节点是异常叶节点
                    % **************** 更新粒球二叉树和队列 ****************
                    [gb1, gb2] = gb_now.self_divide_2_means();
                    gb1.num = tree.node_num + 1;
                    gb2.num = tree.node_num + 2;
                    [tree.node(num_now), tree.node(gb1.num), tree.node(gb2.num)] = gb_now.establish_relationship(gb1, gb2);
                    tree.node_num = tree.node_num + 2;
                    queue = queue.append(gb1.num);
                    queue = queue.append(gb2.num);
                    % (1) 将当前粒球 (gb_now) 分裂为 2 个子粒球 (gb1, gb2)
                    % (2)-(3) 设置当前粒球 (gb_now) 的 2 个子粒球 (gb1, gb2) 的编号
                    % (4) 建立 gb_now 与 gb1 与 gb2 之间的父子关系并在粒球二叉树 (tree) 中更新对应的节点
                    % (5) 更新粒球二叉树中节点的数量 (tree.node_num)
                    % (6)-(7) 将新生成的粒球的编号 (gb1.num, gb2.num) 纳入队列 queue 中
                    % **************** 更新粒球二叉树和队列 ****************
                else % 如果当前叶节点不是异常节点
                    leaf_queue = leaf_queue.append(num_now);
                    % (1) 将当前节点的编号纳入叶节点队列 leaf_queue
                end
            end
            leaf_queue.element = leaf_queue.element(end:-1:1);
            % (1) 将叶节点队列存储的元素 (leaf_queue.element) 换向 (从后向前索引)
            %     其目的是加快处理的速度.
        end

        function tree = update_tree(tree, gamma_range, lambda)
            % 更新粒球二叉树中每个粒球的质量, 初始化最佳质量和最佳子粒球组合
            % param: tree 粒球二叉树
            % param: gamma_range 粒度水平取值范围
            % param: lambda 每个粒球的惩罚系数
            
            % ******************* 初始化检测异常 *******************
            if size(gamma_range, 1) == 0
                fprintf('警告: 粒度水平取值范围为空')
            else
                gamma = gamma_range(1) + eps;
            end
            % (1) 如果 gamma_range 是空数组, 那么输出警告
            % (2) 如果 gamma_range 不是空数组, 那么粒度水平 (gamma) 是粒度水平
            %     范围的最小值 (gamma_range(1)) + 极小的值
            % ******************* 初始化检测异常 *******************
            
            % ******************* 更新粒球二叉树 *******************
            radius_ave = [tree.node(1:tree.node_num).radius_ave];
            coverage   = [tree.node(1:tree.node_num).coverage];
            specificity = 1 ./ (1 + gamma * radius_ave);
            quality     = specificity .* coverage - lambda;
            for item = 1:tree.node_num
                tree.node(item).granularity_level = gamma;
                tree.node(item).specificity       = specificity(item);
                tree.node(item).quality           = quality(item);
                tree.node(item).best_quality      = quality(item);
                tree.node(item).best_combination  = item;
            end
            % (1)-(2)  radius_ave 和 coverage 分别存储每个节点的半径和覆盖率
            % (3)      计算每个节点的特异性
            % (4)      计算每个节点的惩罚质量
            % (5)-(11) 为每个节点的粒度水平, 特异性, 受惩罚质量赋值, 初始化最佳质量和最佳子粒球组合    
            % ******************* 更新粒球二叉树 *******************
        end

        function tree = calculate_gbs_best_combination(tree, leaf_queue)
            % 根据粒球二叉树和叶节点队列计算粒球二叉树根节点的最佳粒球组合
            % param: tree 粒球二叉树
            % param: leaf_queue 粒球二叉树叶节点队列
            % out: gbs 粒球二叉树根节点的最佳粒球组合

            % ************ 计算根节点的最佳粒球组合 ************
            while leaf_queue.now_element_num > 1
                % 当叶节点队列中元素大于 1, 说明粒球二叉树根节点的最佳粒球组合还没有被更新
                num_now = leaf_queue.get();
                leaf_queue = leaf_queue.delete();
                gb_now = tree.node(num_now);
                if ~gb_now.is_leaf
                    continue;
                end
                if mod(num_now, 2) == 0
                    num_brother = num_now + 1;
                else
                    num_brother = num_now - 1;
                end
                gb_brother = tree.node(num_brother);
                % (1)      获取叶节点队列 (leaf_queue) 队首存储的粒球编号 (num_now)
                % (2)      将叶节点队列 (leaf_queue) 的队首存储的编号删除
                % (3)      根据当前粒球编号在粒球二叉树中获取相应的粒球
                % (4)-(6)  当前粒球如果不是叶节点, 说明当前粒球与其兄弟节点已经更新了
                %          它们的父节点, 此时进行下一轮循环, 无需执行后续代码
                % (7)-(8)  如果当前粒球的编号 (num_now) 是偶数, 那么它的兄弟节点编号 (num_brother) 是它的编号 + 1
                % (9)-(11) 如果当前粒球的编号 (num_now) 是奇数, 那么它的兄弟节点编号 (num_brother) 是它的编号 - 1
                % (12)     获取兄弟节点编号 (num_brother) 在粒球二叉树 (tree) 指向的粒球 (gb_brother)
                % 注意: 粒球二叉树根节点编号为 1, 每一个偶数 x 和 x + 1 指向一对兄弟节点

                if ~gb_brother.is_leaf
                    leaf_queue = leaf_queue.append(num_now);
                    continue;
                else
                    num_father = gb_now.num_father;
                    gb_father = tree.node(num_father);
                    tree.node(num_now).is_leaf = false;
                    tree.node(num_brother).is_leaf = false;
                    tree.node(num_father).is_leaf = true;
                    leaf_queue = leaf_queue.append(num_father);
                end
                % (1) 如果兄弟节点 (gb_brother) 不是叶节点
                %     那么无法更新父节点, 因此将当前节点编号重新纳入叶节点队列
                %     不执行后续代码, 进行下一轮循环
                % (2) 如果兄弟节点 (gb_brother) 是叶节点
                %     (a) 获取父节点编号 (num_father)
                %     (b) 获取父节点编号指向的粒球 (gb_father)
                %     (c) 将当前节点 (tree.node(num_now)) 设置为非叶节点
                %     (d) 将兄弟节点 (tree.node(num_brother)) 设置为非叶节点
                %     (e) 将父节点 (tree.node(num_father)) 设置为叶节点
                %     (f) 将父节点编号 (num_father) 纳入叶节点队列

                if gb_father.quality <= gb_now.best_quality + gb_brother.best_quality
                    tree.node(num_father).best_quality = gb_now.best_quality + gb_brother.best_quality;
                    tree.node(num_father).best_combination = cat(1, gb_now.best_combination, gb_brother.best_combination);
                end
                % 如果当前父节点的质量 (gb_father.quality) 小于等于 孩子节点的最佳质量之和 (gb_now.best_quality + gb_brother.best_quality)
                % (1) 更新父节点的最佳质量 (tree.node(num_father).best_quality)
                %     为孩子节点的最佳质量之和 (gb_now.best_quality + gb_brother.best_quality)
                % (2) 更新父节点的最佳子粒球组合 (tree.node(num_father).best_combination)
                %     为孩子节点的最佳子粒球组合的并集 (cat(1, gb_now.best_combination, gb_brother.best_combination))
            end
        end

        function out_gbs = detect_abnormal_gbs(tree, lambda)
            % 输入一些粒球, 检测出相对异常的粒球并分裂它们, 输出无异常的粒球
            % param: in_gbs 输入的粒球 (n*1, GranularBall)
            % out: out_gbs 输出的粒球 (m*1, GranularBall)

            % ************************** 初始化 **************************
            gbs_index = tree.node(1).best_combination;
            value.radius_max = [tree.node(gbs_index).radius_max]';
            value.radius_ave = [tree.node(gbs_index).radius_ave]';
            value.instance_num = [tree.node(gbs_index).instance_num]';
            [value.radius_max_std, mean_value.radius_max_mean] = std(value.radius_max);
            [value.radius_ave_std, mean_value.radius_ave_mean] = std(value.radius_ave);
            [value.instance_num_std, mean_value.instance_num_mean] = std(value.instance_num);
            leaf_queue = Queue(tree.node(1).instance_num, "GranularBall");
            queue = Queue(2 * tree.node(1).instance_num, "uint32");
            % (1) gbs_index 获取根节点的最佳子粒球组合
            % (2)-(4) 获取每个粒球的最大半径、平均半径和样本个数
            % (5)-(7) 计算叶节点最大半径、平均半径以及样本数量的标准差和均值
            % (8) 初始化存储异常节点编号的顺序循环队列 queue
            % (9) 初始化存储生成的粒球的顺序循环队列 leaf_queue
            % ************************** 初始化 **************************

            % ************************ 计算异常节点 ************************
            abnormal_max_index = value.radius_max > (mean_value.radius_max_mean + value.radius_max_std);
            abnormal_ave_index = (value.radius_ave > (mean_value.radius_ave_mean + value.radius_ave_std)) & ...
                (value.instance_num < (mean_value.instance_num_mean - value.instance_num_std));
            abnormal_index = abnormal_ave_index | abnormal_max_index;
            queue.element = gbs_index(abnormal_index);
            queue.now_element_num = sum(abnormal_index);
            queue.rear = queue.now_element_num + 1;
            leaf_queue.element = tree.node(gbs_index(~abnormal_index));
            leaf_queue.now_element_num = sum(~abnormal_index);
            leaf_queue.rear = leaf_queue.now_element_num + 1;
            % (1) - (3) 获取异常叶节点的索引 abnormal_index
            % (4) - (6) 初始化存储异常节点编号的队列 queue
            % (7) - (9) 初始化存储正常粒球的队列 leaf_queue

            while ~queue.is_empty() % 持续检测直到所有粒球都被检测了
                num_now = queue.get();
                queue = queue.delete();
                gb_now = tree.node(num_now);
                % (1) 从队列 (queue) 中获取队首的粒球的编号
                % (2) 删除队列的队首
                % (3) 获取编号对应的粒球

                if gb_now.num_left_child ~= -1
                    % 如果当前粒球有左右孩子节点
                    gb1 = tree.node(gb_now.num_left_child);
                    gb2 = tree.node(gb_now.num_right_child);
                    % (1)-(2) gb1 和 gb2 分别 gb_now 的左右孩子节点
                else
                    % 如果当前粒球没有左右孩子节点
                    [gb1, gb2] = gb_now.self_divide_2_means();
                    gb1.num = tree.node_num + 1;    gb2.num = tree.node_num + 2;
                    gb1.granularity_level = gb_now.granularity_level;   gb2.granularity_level = gb_now.granularity_level;
                    gb1.specificity = 1 / (1 + gb1.radius_ave * gb1.granularity_level); gb2.specificity = 1 / (1 + gb2.radius_ave * gb2.granularity_level);
                    gb1.quality = gb1.coverage * gb1.specificity - lambda;  gb2.quality = gb2.coverage * gb2.specificity - lambda;
                    [tree.node(num_now), tree.node(gb1.num), tree.node(gb2.num)] = gb_now.establish_relationship(gb1, gb2);
                    tree.node_num = tree.node_num + 2;
                    % (1) 分裂当前粒球得到 gb1 和 gb2
                    % (2) 为 gb1 和 gb2 设置节点编号
                    % (3) 为 gb2 和 gb2 设置粒度水平
                    % (4) 计算 gb1 和 gb2 的特异性
                    % (5) 计算 gb1 和 gb2 的质量
                    % (6) 为 gb_now, gb1 和 gb2 建立关系, 并将其赋值至粒球二叉树
                    % (7) 粒球二叉树节点数量 + 2
                end
                if gb1.radius_max > mean_value.radius_max_mean + value.radius_max_std || ...
                        (gb1.radius_ave > mean_value.radius_ave_mean + value.radius_ave_std && ...
                        gb1.instance_num < mean_value.instance_num_mean - value.instance_num_std)
                    queue = queue.append(gb1.num);
                else
                    leaf_queue = leaf_queue.append(gb1);
                end
                if gb2.radius_max > mean_value.radius_max_mean + value.radius_max_std || ...
                        (gb2.radius_ave > mean_value.radius_ave_mean + value.radius_ave_std && ...
                        gb2.instance_num < mean_value.instance_num_mean - value.instance_num_std)
                    queue = queue.append(gb2.num);
                else
                    leaf_queue = leaf_queue.append(gb2);
                end
                % (1) 如果 gb1 是异常的节点, 那么将 gb1 的编号纳入异常粒球队列 queue
                % (2) 如果 gb1 不是异常节点, 那么将 gb1 纳入正常粒球队列 leaf_queue
                % (3) 如果 gb2 是异常的节点, 那么将 gb2 的编号纳入异常粒球队列 queue
                % (4) 如果 gb2 不是异常节点, 那么将 gb2 纳入正常粒球队列 leaf_queue
            end
            out_gbs = leaf_queue.element(1:leaf_queue.rear-1);
            % 从 leaf_queue 将输出的粒球构成的数组提取出来
        end
    end
end

