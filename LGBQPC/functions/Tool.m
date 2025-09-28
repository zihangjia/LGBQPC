classdef Tool
    %TOOL 工具类
    % 提供以下静态方法:
    %
    % method 1: solution = solve_3_order_polynomial_inequality(in_param)
    %   解最高三次多项式的不等式 f(x) = ax^3 + bx^2 + cx + d >= 0
    %   输入:  in_param = [a, b, c, d]
    %   输出:  解区间矩阵 solution (k x 2), 每行表示一个区间
    %
    % method 2: out_intervals = interval_intersection(in_intervals_1, in_intervals_2)
    %   计算两个区间集合的交集
    %   输入:  in_intervals_1, in_intervals_2 (n x 2 矩阵, 每行一个区间)
    %   输出:  out_intervals (m x 2 矩阵, 表示交集区间集合)
    %
    % method 3: [is_not_spherical, p_value, chi_square_statistic, degree_freedom] = ...
    %            bartlett_sphericity_test(data_matrix, significance_level)
    %   进行 Bartlett 球形性检验
    %   输入:  data_matrix (n*m 样本矩阵), significance_level 显著性水平
    %   输出:  is_not_spherical      是否拒绝原假设 (1=非球形, 0=球形)
    %          p_value              卡方检验的 p 值
    %          chi_square_statistic 卡方统计量
    %          degree_freedom       自由度 p*(p-1)/2

    methods(Static)
        function solution = solve_3_order_polynomial_inequality(in_param)
            %SOLVE_3_ORDER_POLYNOMIAL_INEQUALITY
            % 解一元三次及以下多项式不等式 f(x) >= 0
            %
            % 输入:
            %   in_param = [a, b, c, d] 表示 f(x) = a*x^3 + b*x^2 + c*x + d
            %
            % 输出:
            %   solution 为 k x 2 的矩阵, 每行表示解集中的一个区间
            %   例如: [-inf, r1; r2, r3] 表示解为 (-inf, r1] ∪ [r2, r3]

            % *************************** 初始化 ***************************
            tolerance = 1e-10;
            a = in_param(1); b = in_param(2); c = in_param(3); d = in_param(4);
            % (1) 初始化误差容忍度 tolerance
            % (2) 初始化 3 次不等式的系数
            % *************************** 初始化 ***************************

            % *********************** 判断不等式次数 **********************
            if abs(a) > tolerance
                polynomial_degree = 3; root = roots([a, b, c, d]);
            elseif abs(b) > tolerance
                polynomial_degree = 2; root = roots([b c d]);
            elseif abs(c) > tolerance
                polynomial_degree = 1; root = -d/c;
            else
                polynomial_degree = 0;
            end
            % *********************** 判断不等式次数 **********************


            % ***************** 根据多项式的次数求解不等式 *****************
            switch polynomial_degree
                case 0 % 常函数 f(x) = d
                    if d >= 0 % f(x) >= 0 恒成立
                        solution = [-inf, inf];
                    else % f(x) >= 0 恒不成立
                        solution = [];
                    end
                case 1 % 一次函数 f(x) = c*x + d
                    if c > 0 % 增函数, x >= root 时, f(x) >= 0
                        solution = [root, inf];
                    else % 减函数, x <= root 时, f(x) <= 0
                        solution = [-inf, root];
                    end
                case 2 % 二次函数 f(x) = b*x^2 + c*x + d
                    root = unique(real(root(abs(imag(root))<tolerance))); 
                    % 如果根的虚部的绝对值超过容忍度将会被剔除, 仅保留不同根的实部并对根进行排序
                    root_num = numel(root); % 实根的个数
                    switch root_num
                        case 0 % 如果没有实根
                            if b > 0 % 二次函数开口向上，始终 > 0
                                solution = [-inf, inf];
                            else % 二次函数开口向下，始终 < 0
                                solution = [];
                            end
                        case 1
                            if b > 0 % 二次函数开口向上
                                solution = [-inf, inf];
                            else % 开口向下, 解集为两根之间的区域
                                solution = [root(1), root(1)]; % 开口向下 → 夹在根之间
                            end
                        case 2 % 如果有两个实根, 包含重根
                            if b > 0 % 二次函数开口向上
                                solution = [-inf, root(1); root(2), inf];
                            else % 开口向下, 解集为两根之间的区域
                                solution = [root(1), root(2)]; % 开口向下 → 夹在根之间
                            end
                    end
                case 3 % 三次函数 f(x) = a*x^3 + b*x^2 + c*x + d
                    root = unique(real(root(abs(imag(root))<tolerance))); 
                    % 如果根的虚部的绝对值超过容忍度将会被剔除, 仅保留不同根的实部并对根进行排序
                    root_num = numel(root); % 实根的个数
                    switch root_num
                        case 0 % 没有实根
                            if a > 0 % a > 0 解集是全空间
                                solution = [-inf, inf];
                            else % a < 0 解集是空集
                                solution = [];
                            end
                        case 1 % 1 个实根
                            if a > 0 % a > 0 解集是根到正无穷
                                solution = [root(1), inf];
                            else % a > 0 解集是负无穷到根
                                solution = [-inf, root(1)];
                            end
                        case 2 % 2 个实根
                            root_mid = (root(1) + root(2)) / 2; % 两根的中点
                            if a > 0 % f(x)->inf, x->inf
                                if a * power(root_mid, 3) + b * power(root_mid, 2) + c * root_mid + d > 0
                                    solution = [root(1), inf];
                                else
                                    solution = [root(1), root(1); root(2), inf];
                                end
                            else % f(x)->-inf, x->inf
                                if a * power(root_mid, 3) + b * power(root_mid, 2) + c * root_mid + d > 0
                                    solution = [-inf, root(2)];
                                else
                                    solution = [-inf, root(1); root(2), root(2)];
                                end
                            end
                        case 3
                            if a > 0 % f(x)->inf, x->inf
                                solution = [root(1), root(2); root(3), inf];
                            else % f(x)->-inf, x->inf
                                solution = [-inf, root(1); root(2), root(3)];
                            end
                    end
            end
            % ***************** 根据多项式的次数求解不等式 *****************
        end

        function out_intervals = interval_intersection(in_intervals_1, in_intervals_2)
            %INTERVAL_INTERSECTION 计算两个区间集合的交集
            %
            % 输入:
            %   in_intervals_1 - n1 x 2 矩阵, 每行一个区间 [l, r]
            %   in_intervals_2 - n2 x 2 矩阵
            %
            % 输出:
            %   out_intervals  - m x 2 矩阵, 表示交集区间集合

            % ************************** 输入检查 **************************
            if isempty(in_intervals_1) || isempty(in_intervals_2)
                out_intervals = zeros(0,2);
                return;
            end
            in_intervals_1 = sortrows(in_intervals_1, 1);
            in_intervals_2 = sortrows(in_intervals_2, 1);
            % (1)-(4) 如果区间集合 in_intervals_1 或 in_intervals_2 为空, 则返回空集
            % (5)-(6) 对区间集合按照左端点升序排列
            % ************************** 输入检查 **************************

            % 确保输入已排序 (按区间左端点升序)


            % ************************ 双指针遍历 ************************
            item_1 = 1; item_2 = 1;
            intervals_1_num = size(in_intervals_1, 1);
            intervals_2_num = size(in_intervals_2, 1);
            results = zeros(intervals_1_num + intervals_2_num, 2);
            count = 0;
            % (1) 初始化指针 item_1 和 item_2
            % (2)-(3) 计算区间集合 intervals_1 和 intervals_2 中区间的个数
            % (4)-(5) 初始化区间集合交的结果和非空交集的数量
            while item_1 <= intervals_1_num && item_2 <= intervals_2_num
                left  = max(in_intervals_1(item_1, 1), in_intervals_2(item_2, 1));
                right = min(in_intervals_1(item_1, 2), in_intervals_2(item_2, 2));
                if left <= right
                    count = count + 1;
                    results(count,:) = [left, right];
                end
                if in_intervals_1(item_1,2) < in_intervals_2(item_2,2)
                    item_1 = item_1 + 1;
                else
                    item_2 = item_2 + 1;
                end
                % (1)-(2)  计算两个区间交集的左端点和右端点
                % (3)-(6)  如果左端点小于右端点, 说明有交集. 非空交集的数量 + 1
                %          且第 count 个非空交集的结果是当前的交集
                % (7)-(11) 如果区间 in_intervals_1(item_1, :) 的右端点小, 说明
                %          它和区间集合 in_intervals_2 中 item_2 以后的区间无交集
                %          所以指针 item_1 指向下一个区间
            end
            % ************************ 双指针遍历 ************************

            % ************************ 输出结果 ************************
            out_intervals = results(1:count,:);
            % 截断无效部分
            % ************************ 输出结果 ************************
        end


        function [is_not_spherical, p_value, chi_square_statistic, degree_freedom] = ...
                bartlett_sphericity_test(data_matrix, significance_level)
            %BARTLETTSPHERICITYTEST Bartlett 球形性检验
            %
            % 输入:
            %   data_matrix        - n x m 数值矩阵 (n 个样本, m 个变量)
            %   significance_level - 显著性水平 (可选, 默认 0.01)
            %
            % 输出:
            %   is_not_spherical       - 是否拒绝原假设 (1=非球形, 0=球形)
            %   p_value                - p 值 (若无法计算则返回 NaN)
            %   chi_square_statistic   - 卡方统计量 (若无法计算则 NaN)
            %   degree_freedom         - 自由度 = p*(p-1)/2

            % ********************** 默认参数设置 **********************
            if nargin < 2 || isempty(significance_level)
                significance_level = 0.01;
            end
            % (1)-(3) 如果未设置显著性水平, 则设置显著性水平为 0.01
            % ********************** 默认参数设置 **********************

            % ************************ 输入检查 ************************
            if ~ismatrix(data_matrix) || ~isnumeric(data_matrix)
                error('data_matrix 必须是数值型二维矩阵。');
            elseif any(isnan(data_matrix(:)))
                error('data_matrix 包含 NaN，请先处理缺失值。');
            end
            [instance_num, feature_num] = size(data_matrix);
            if feature_num < 2
                error('变量个数 m 必须 >= 2 才能进行 Bartlett 检验。');
            end
            % (1)-(2) 检查 data_matrix 是否是数值型二维矩阵
            % (3)-(5) 检查 data_matrix 是否含有缺失值, matrix(:) 是列向量, isnan()
            %         检查是否有缺失值, any() 检查输入的向量或矩阵的第 1 列是否有非零值或 true
            % (6)     instance_num 和 feature_num 存储样本数量和特征数量
            % (7)-(9) 检查变量的个数是否大于 2
            % ************************ 输入检查 ************************


            % ******************** 判断能否进行统计检验 ********************
            correlation_matrix = corr(data_matrix);
            correlation_matrix = (correlation_matrix + correlation_matrix') / 2;
            eig_values = eig(correlation_matrix);
            eig_values_min = min(eig_values);
            degree_freedom = feature_num * (feature_num - 1) / 2;
            if eig_values_min <= 0
                is_not_spherical = true;
                p_value = NaN;
                chi_square_statistic = NaN;
                warning('相关矩阵含非正特征值 (minEig = %.3e)。', eig_values_min);
                return;
            end
            % (1)-(2)  计算相关矩阵并强制让相关矩阵对称
            % (3)      eig_values 存储相关矩阵的特征值
            % (4)      eig_values_min 存储相关矩阵的特征值的最小值
            % (5)      计算自由度
            % (6)-(12) 如果相关矩阵的特征值的最小值小于等于 0, 说明相关矩阵是奇
            %          异矩阵无法进行统计检验, 此时数据集一定非球形, p 值和卡方
            %          统计量设置为 NaN. 最后, 输出警告并返回
            % ******************** 判断能否进行统计检验 ********************

            % ************************ 计算统计量 ************************
            logDetR = sum(log(eig_values));
            c = (instance_num - 1 - (2*feature_num + 5)/6);
            if c <= 0
                warning('样本数太小，修正系数 c = %.4f ≤ 0，近似可能不可靠。', c);
            end
            chi_square_statistic = - c * logDetR;
            % (1) 用特征值计算相关系数矩阵的行列式的对数
            % (2) 计算统计量的修正系数
            % (3)-(5) 如果统计量的修正系数太小, 发出警告
            % (6) 计算卡方统计量
            % ************************ 计算统计量 ************************

            % ********************* 给出统计检验结果 *********************
            if ~isreal(chi_square_statistic) || ~isfinite(chi_square_statistic)
                is_not_spherical = true;
                p_value = NaN;
                chi_square_statistic = NaN;
                warning('卡方统计量计算异常。');
                return;
            end
            p_value = 1 - chi2cdf(chi_square_statistic, degree_freedom);
            is_not_spherical = (p_value < significance_level);
            % (1)-(7) 如果卡方统计量不是实数或者不是有限数, 则认为数据集非球形
            %         p 值和卡方统计量是 NaN, 并输出警告
            % (8) 计算 p 值
            % (9) 如果 p 值小于显著性水平, 则认为数据集非球形, 否则数据集是球形
            % ********************* 给出统计检验结果 *********************
        end
    end
end
