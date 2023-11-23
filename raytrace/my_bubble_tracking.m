function bubble_tracking_out = my_bubble_tracking(bubble,varargin)
% 本函数的作用是为气泡图像处理算法提供二维、三维跟踪功能，使用方式如下：
% 
% 输入：
%     bubble：气泡信息  提取：（拟合所用边界/质心/中心/尺寸）
%     varargin：可变长度输入列表
% 
% 输出：
%     bubble_tracking_out：气泡追踪结果 （该气泡开始存在时间/消失时间/存在期间信息{存在1{拟合所用边界{[x,x]...}/质心{[x,x]}/中心{[x,x]}/尺寸{[x,x,x]}}/存在2/...}）
% 
% 调用说明：
%     bubble_tracking_out = my_bubble_tracking(bubble,'first')
%                  第一时刻图的输入方式
%     bubble_tracking_out = my_bubble_tracking(bubble,'others',bubble_tracking_data)
%                  其他时刻图的输入方式
%     bubble_tracking_out = my_bubble_tracking(bubble,'first'/'others',varargin)
%                  可识别关键词：  'maxlength_th'（最大质心差） '3_d'（3维气泡追踪）
% 
% 版本号VOL1.0，编写于2021年7月8日，作者：WG-Chen
% 重要：要注意识别过程如果找到了同一个气泡的两个（或以上）下一时刻的气泡，的判别方法：1对比质心差、2对比尺寸信息，综合对比

%% default value
flag = 2 ;          % 默认二维重建
maxlength_th = 5;   % 默认最大质心差为5
method = 1;         % 默认输入第一时刻图像
u_v_flag = 0;       % 默认不计算气泡速度

%% methods
if numel(varargin) > 0
    tf = {'first','others','maxlength_th','3_d','fps'};
    for i = 1 : numel(tf)
        for j = 1 : numel(varargin)
            tf_1 = strcmp(tf{i},varargin{j});
            if tf_1 == 1
                switch i
                    case 2
                        method = 2;
                        bubble_tracking_data = varargin{j+1};
                    case 3
                        maxlength_th = varargin{j+1};
                    case 4
                        flag = 3;
                    case 5
                        u_v_flag = 1;
                        fps = varargin{j+1};
                end
            end
        end
    end
end

%% main

if flag == 2 
    if method == 1
        bubble_num = 0;
        bubble_tracking_out = {};
        [num,~] = size(bubble);
        for i = 1 : num
            if bubble{i,2} > 1
                [b,~]=size(bubble{i,6});
                for j = 1 : b
                    if bubble{i,6}{j,1} ~= 0
                        bubble_num = bubble_num + 1;
                        bubble_tracking_out{bubble_num,1} = 1;
                        bubble_tracking_out{bubble_num,2} = 1;
                        bubble_tracking_out{bubble_num,3}{1,1}{1,1} = bubble{i,7}{j,1};
                        bubble_tracking_out{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7}{j,1});
                        bubble_tracking_out{bubble_num,3}{1,1}{1,3} = [bubble{i,6}{j,1}(1,3),bubble{i,6}{j,1}(1,4)];
                        bubble_tracking_out{bubble_num,3}{1,1}{1,4} = [bubble{i,6}{j,1}(1,1),bubble{i,6}{j,1}(1,2),bubble{i,6}{j,1}(1,5)];
                    end
                end
            elseif bubble{i,2} == 1 || bubble{i,2} == 0
                if bubble{i,6} ~= 0
                    bubble_num = bubble_num + 1;
                    bubble_tracking_out{bubble_num,1} = 1;
                    bubble_tracking_out{bubble_num,2} = 1;
                    bubble_tracking_out{bubble_num,3}{1,1}{1,1} = bubble{i,7};
                    bubble_tracking_out{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7});
                    bubble_tracking_out{bubble_num,3}{1,1}{1,3} = [bubble{i,6}(1,3),bubble{i,6}(1,4)];
                    bubble_tracking_out{bubble_num,3}{1,1}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,2),bubble{i,6}(1,5)];
                end
            else
                if bubble{i,6} ~= 0
                    bubble_num = bubble_num + 1;
                    bubble_tracking_out{bubble_num,1} = 1;
                    bubble_tracking_out{bubble_num,2} = 1;
                    bubble_tracking_out{bubble_num,3}{1,1}{1,1} = bubble{i,7};
                    bubble_tracking_out{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7});
                    bubble_tracking_out{bubble_num,3}{1,1}{1,3} = [bubble{i,6}(2,1),bubble{i,6}(2,2)];
                    bubble_tracking_out{bubble_num,3}{1,1}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,1),0];
                end
            end
        end
        for kkk = 1 : bubble_num
            bubble_tracking_out{kkk,4} = 0;
        end
    else
        [num,~] = size(bubble);
        [bubble_num,~] = size(bubble_tracking_data);
        bubble_num_raw = bubble_num;
        monitor = zeros(bubble_num_raw,3);  % 1 是否有下一个气泡  / [2,3] 下一个气泡在bubble上的位置
        for i = 1 : num
            if bubble{i,2} > 1
                [b,~]=size(bubble{i,6});
                for j = 1 : b
                    if bubble{i,6}{j,1} ~= 0
                        for n = 1 : bubble_num_raw
                            leng = norm([bubble{i,6}{j,1}(1,3),bubble{i,6}{j,1}(1,4)]-bubble_tracking_data{n,3}{1,bubble_tracking_data{n,1}}{1,3});
                            if leng <= maxlength_th & monitor(n,1) == 0 & bubble_tracking_data{n,4} ~= 1
                                monitor(n,1) = 1;monitor(n,2) = i;monitor(n,3) = j; % 监视器
                                bubble{i,6}{j,1}(1,6) = 1;% 监视器
                                bubble_tracking_data{n,1} = bubble_tracking_data{n,1}+1;
                                haha = bubble_tracking_data{n,1};
                                bubble_tracking_data{n,3}{1,haha}{1,1} = bubble{i,7}{j,1};
                                bubble_tracking_data{n,3}{1,haha}{1,2} = my_centeroid(bubble{i,7}{j,1});
                                bubble_tracking_data{n,3}{1,haha}{1,3} = [bubble{i,6}{j,1}(1,3),bubble{i,6}{j,1}(1,4)];
                                bubble_tracking_data{n,3}{1,haha}{1,4} = [bubble{i,6}{j,1}(1,1),bubble{i,6}{j,1}(1,2),bubble{i,6}{j,1}(1,5)];
                                if u_v_flag
                                    bubble_tracking_data{n,4}(bubble_tracking_data{n,1}-1,1) = leng*fps;  % unit=(pxs/s)
                                end
                                break
                            end
                        end
                    end
                end
            elseif bubble{i,2} == 0 || bubble{i,2} == 1
                if bubble{i,6} ~= 0
                    for n = 1 : bubble_num_raw
                        leng = norm([bubble{i,6}(1,3),bubble{i,6}(1,4)]-bubble_tracking_data{n,3}{1,bubble_tracking_data{n,1}}{1,3});
                        if leng <= maxlength_th & monitor(n,1) == 0 & bubble_tracking_data{n,4} ~= 1
                            monitor(n,1) = 1;monitor(n,2) = i;monitor(n,3) = 1; % 监视器
                            bubble{i,6}(1,6) = 1;% 监视器
                            bubble_tracking_data{n,1} = bubble_tracking_data{n,1}+1;
                            haha = bubble_tracking_data{n,1};
                            bubble_tracking_data{n,3}{1,haha}{1,1} = bubble{i,7};
                            bubble_tracking_data{n,3}{1,haha}{1,2} = my_centeroid(bubble{i,7});
                            bubble_tracking_data{n,3}{1,haha}{1,3} = [bubble{i,6}(1,3),bubble{i,6}(1,4)];
                            bubble_tracking_data{n,3}{1,haha}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,2),bubble{i,6}(1,5)];
                            if u_v_flag
                                bubble_tracking_data{n,4}(bubble_tracking_data{n,1}-1,1) = leng*fps;  % unit=(pxs/s)
                            end
                            break
                        end
                    end
                end
            else
                if bubble{i,6} ~= 0
                    for n = 1 : bubble_num_raw
                        leng = norm([bubble{i,6}(2,1),bubble{i,6}(2,2)]-bubble_tracking_data{n,3}{1,bubble_tracking_data{n,1}}{1,3});
                        if leng <= maxlength_th & monitor(n,1) == 0 & bubble_tracking_data{n,4} ~= 1
                            monitor(n,1) = 1;monitor(n,2) = i;monitor(n,3) = 1; % 监视器
                            bubble{i,6}(1,6) = 1;% 监视器
                            bubble_tracking_data{n,1} = bubble_tracking_data{n,1}+1;
                            haha = bubble_tracking_data{n,1};
                            bubble_tracking_data{n,3}{1,haha}{1,1} = bubble{i,7};
                            bubble_tracking_data{n,3}{1,haha}{1,2} = my_centeroid(bubble{i,7});
                            bubble_tracking_data{n,3}{1,haha}{1,3} = [bubble{i,6}(2,1),bubble{i,6}(2,2)];
                            bubble_tracking_data{n,3}{1,haha}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,1),0];
                            if u_v_flag
                                bubble_tracking_data{n,4}(bubble_tracking_data{n,1}-1,1) = leng*fps;  % unit=(pxs/s)
                            end
                            break
                        end
                    end
                end
            end
        end
        % 处理原来没有的，即新增/原来有的，没有了，即消失（离焦、离开边界）
        for i = 1 : bubble_num_raw
            if monitor(i,1) == 0
                bubble_tracking_data{i,2} = bubble_tracking_data{i,1};
                bubble_tracking_data{i,4} = 1;
            end
        end
        
        for i = 1 : num
            if bubble{i,2} > 1
                [b,~]=size(bubble{i,6});
                for j = 1 : b
                    if bubble{i,6}{j,1} ~= 0 & length(bubble{i,6}{j,1}) == 5
                        bubble_num = bubble_num + 1;
                        bubble_tracking_data{bubble_num,1} = 1;
                        bubble_tracking_data{bubble_num,2} = 1;
                        bubble_tracking_data{bubble_num,3}{1,1}{1,1} = bubble{i,7}{j,1};
                        bubble_tracking_data{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7}{j,1});
                        bubble_tracking_data{bubble_num,3}{1,1}{1,3} = [bubble{i,6}{j,1}(1,3),bubble{i,6}{j,1}(1,4)];
                        bubble_tracking_data{bubble_num,3}{1,1}{1,4} = [bubble{i,6}{j,1}(1,1),bubble{i,6}{j,1}(1,2),bubble{i,6}{j,1}(1,5)];
                    end
                end
            elseif bubble{i,2} == 1 || bubble{i,2} == 0
                if bubble{i,6} ~= 0 & length(bubble{i,6}) == 5
                    bubble_num = bubble_num + 1;
                    bubble_tracking_data{bubble_num,1} = 1;
                    bubble_tracking_data{bubble_num,2} = 1;
                    bubble_tracking_data{bubble_num,3}{1,1}{1,1} = bubble{i,7};
                    bubble_tracking_data{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7});
                    bubble_tracking_data{bubble_num,3}{1,1}{1,3} = [bubble{i,6}(1,3),bubble{i,6}(1,4)];
                    bubble_tracking_data{bubble_num,3}{1,1}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,2),bubble{i,6}(1,5)];
                end
            else
                if bubble{i,6} ~= 0 & length(bubble{i,6}) == 2
                    bubble_num = bubble_num + 1;
                    bubble_tracking_data{bubble_num,1} = 1;
                    bubble_tracking_data{bubble_num,2} = 1;
                    bubble_tracking_data{bubble_num,3}{1,1}{1,1} = bubble{i,7};
                    bubble_tracking_data{bubble_num,3}{1,1}{1,2} = my_centeroid(bubble{i,7});
                    bubble_tracking_data{bubble_num,3}{1,1}{1,3} = [bubble{i,6}(2,1),bubble{i,6}(2,2)];
                    bubble_tracking_data{bubble_num,3}{1,1}{1,4} = [bubble{i,6}(1,1),bubble{i,6}(1,1),0];
                end
            end
        end
        for kkk = bubble_num_raw+1 : bubble_num
            bubble_tracking_data{kkk,4} = 0;
        end
        bubble_tracking_out = bubble_tracking_data;
    end
end



end

function out = my_centeroid(bubble_bw)

    sumx = 0;
    sumy = 0;
    
    for j = 1 : size(bubble_bw,1)
        sumx = sumx + bubble_bw(j,1);
        sumy = sumy + bubble_bw(j,2);
    end
    
    xc = sumx/size(bubble_bw,1);
    yc = sumy/size(bubble_bw,1);
    
    out = {};
    out{1,1}(1,1) = xc;
    out{1,1}(1,2) = yc;

end
