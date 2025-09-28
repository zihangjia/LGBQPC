classdef QueueNum
    % Queue 实现先进先出队列类
    %   此处显示详细说明
    
    properties
        element         % 存储元素的元胞数组
        front           % 队首指标
        rear            % 队尾指标
        now_element_num % 队列中当前元素个数
        max_element_num % 队列中最大元素个数
    end
    
    methods
        function obj = QueueNum(input_max_element_num)
            % 实例化一个队列
            % param: input_max_element_num 队列中所含最大元素个数
            
            obj.element = -1 * ones(input_max_element_num, 1);
            % 队列用 input_max_element_num * 1 的元胞数组存储元素
            obj.max_element_num = input_max_element_num;
            % 初始化队列最大元素个数
            obj.front = 1;
            obj.rear = 1;
            obj.now_element_num = 0;
            % 初始化队列的首 (front) 尾 (rear) 指向 1, 当前元素个数为 0
        end
        
        function output = is_empty(obj)
            % 如果当前队列为空返回 true, 否则返回 false
            
            if obj.now_element_num == 0
                output = true;
            else
                output = false;
            end
        end

        function output = get(obj)
            % 获取队列的队首元素

            % 当队列中元素的数量不为 0 的情况下返回队列中 front 指向的元素 
            % 否则输出队列已空
            if obj.now_element_num ~= 0 
                output = obj.element(obj.front);
            else
                fprintf('队列已空, 无法获取元素\n');
            end
        end

        function obj = append(obj, input_element)
            % 在队列中插入输入的元素
            % param: input_element 输入的元素

            if obj.now_element_num > 0 && obj.front == obj.rear
                % 如果队首 (front) 和队尾 (rear) 指向同一个元素且队列不为空
                % 则顺序循环队列已满
                fprintf('队列已满无法插入！\n');
            else % 否则插入 dataset
                obj.element(obj.rear) = input_element;
                % 将元素放置在队尾 (rear) 指向的位置
                obj.rear = mod(obj.rear + 1, obj.max_element_num + 1);
                % 如果 obj.rear + 1 指向 obj.max_element+1, 则暂设 obj.rear 为 0
                % 否则 obj.rear = obj.rear + 1
                if obj.rear == 0 
                % 因为 MATLAB 中数组从 1 开始, 所以 obj.rear 是 0 的话调整为 1
                    obj.rear = obj.rear + 1;
                end
                obj.now_element_num = obj.now_element_num + 1;
                % 队列中元素个数 + 1
            end
        end

        function obj = delete(obj)
            % 在顺序循环队列中删除队首元素 

            if obj.now_element_num == 0
                % 如果顺序循环队列中已经没有元素, 则什么也不做
                fprintf('队列已空无数据元素出队列！\n');
            else % 否则, 从顺序循环队列中删除队首指标指向的元素
                obj.front = mod(obj.front + 1, obj.max_element_num + 1);
                % 如果 obj.front + 1 指向 obj.max_element+1, 则暂设 obj.front 为 0
                % 否则 obj.front = obj.front + 1
                if obj.front == 0
                % 因为 matlab 中数组从 1 开始, 所以 obj.front 是 0 的话调整为 1
                    obj.front = obj.front + 1;
                end
                obj.now_element_num = obj.now_element_num - 1;
                % 队列中元素个数 - 1
            end
        end
    end
end

