function bubble_out = my_statistics(bubble_in,varargin)
% 本函数的作用是为气泡图像处理算法提供统计分析，将离平均气泡尺寸很远的气泡删除，使用方式如下：
% 本函数是假设气泡尺寸在一次实验中满足对数正态分布，并认为大于上分位点α=0.05的气泡属于小概率事件，认定为识别误差，予以去除； 
% 由于在实际操作中，使用样本均值和方差来估计总体均值和方差，应要求气泡的数量满足N >> 45；所以当不满足时，不予修正，返回原来的气泡数据；
%
% 输入：
%     bubble_in：输入的气泡尺寸信息（1×N的数组）
%     varargin：可变输入列表，用于控制阈值
% 
% 输出：
%     bubble_out：经过校正后的气泡尺寸信息（1×N的数组）
% 
% 调用说明：
%     bubble_out = my_statistics(bubble_in)
%                  默认预处理方式
%     bubble_out = my_statistics(bubble_in,0.1,100)
%                  设定α值为0.1，默认值为0.05;设定气泡数量阈值为100，默认值为60；
% 
% 版本号VOL1.0，编写于2021年6月6日，作者：WG-Chen

%% default value
alpha = 0.05;    %似乎太大了，设置为0.02较为合适
num_th = 45;

%% method
if numel(varargin) > 0
    for i = 1 : numel(varargin)
        if varargin{i} < 1
            alpha = varargin{i};
        else
            num_th = varargin{i};
        end
    end
end


%% calculating

[Num,~] = size(bubble_in);
if Num > num_th
    E_X = mean(bubble_in);
    D_X = var(bubble_in);

    mu = log(E_X)-0.5*log(1+D_X/E_X^2);
    sigma = sqrt(log(1+D_X/E_X^2));

    p = 1 - alpha;

    bub_size_th = logninv(p,mu,sigma);

    num = 0;
    bubble_out = zeros();
    for i = 1 : length(bubble_in)
        if bubble_in(i,1) <= bub_size_th
            num = num + 1;
            bubble_out(num,1) = bubble_in(i,1);
        end
    end
else
    bubble_out = bubble_in;
end


bubble_out(1,2) = mean(bubble_out(:,1));
[bubble_out(1,3),~] = size(bubble_out);
bubble_out(1,4) = max(bubble_out(:,1));

bubble_out(3,2) = mean(bubble_in(:,1));
[bubble_out(3,3),~] = size(bubble_in);
bubble_out(3,4) = max(bubble_in);

end