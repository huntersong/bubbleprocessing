function bubble_out = my_bub_deleting(bubble,varargin)
% 本函数的作用是为气泡图像处理算法提供气泡识别结果进行校正功能，对单气泡进行分类拟合处理，使用方式如下：
% 
% 输入：
%     bubble：气泡边界信息
%     varargin：可变长度输入列表，用于指定阈值
% 
% 输出：
%     bubble_out：修正后的识别到的气泡信息，包括圆度、质心、尺寸、角度等
% 
% 调用说明：
%     bubble_out = my_bub_deleting(bubble)
%                  默认预处理方式，输出修正后的气泡信息
%     bubble_out = my_bub_deleting(bubble,varargin)
%                  可识别关键词：  'maxsize_th'（最大尺寸） 'minsize_tht'（最小尺寸）
%                  'ratio_th' （椭圆度阈值）'repetition_length_th'（重复椭圆圆心距阈值）
%                  'method'（方法：'repetition','abnormal','total'）
% 
% 版本号VOL1.0，编写于2021年6月4日，作者：WG-Chen
% 版本号VOL2.0，编写于2021年6月29日，作者：WG-Chen
% 修改内容：增加一种气泡删除模块，默认需要进行，优先级为最高，即分析拟合出来的椭圆和
%           拟合使用的弧的拟合度需要大于阈值，默认阈值为0.65

%% default value
maxsize_th = 60;
minsize_th = 5;
ratio_th = 0.4;
maxratio_th = 1/ratio_th;
repetition_length_th = 10;
method_flag = 0;
fitting_th = 0.65;

%% methods
if numel(varargin) > 0
    tf = {'maxsize_th','minsize_th','ratio_th','repetition_length_th','method','fitting_th'};
    for i = 1 : numel(tf)
        for j = 1 : numel(varargin)
            tf_1 = strcmp(tf{i},varargin{j});
            if tf_1 == 1
                switch i
                    case 1
                        maxsize_th = varargin{j+1};
                    case 2
                        minsize_th = varargin{j+1};
                    case 3
                        ratio_th = varargin{j+1};
                    case 4
                        repetition_length_th = varargin{j+1};
                    case 5
                        tf1 = {'abnormal','repetition','total'};
                        tf2 = strcmp(varargin{j+1},tf1);
                        tf3 = find(tf2==1);
                        switch tf3
                            case 1 
                                method_flag = 1;
                            case 2
                                method_flag = 2;
                        end
                    case 6
                        fitting_th = varargin{j+1};
                end
            end
        end
    end
end

%% 计算部
if method_flag == 0
    bubble = fitting_del(bubble,fitting_th);
    bubble = abnormal_ellipse(bubble,maxsize_th,minsize_th,ratio_th,maxratio_th);
    bubble_out = repetitive_ellipse(bubble,repetition_length_th);
elseif method_flag == 1
    bubble_out = abnormal_ellipse(bubble,repetition_length_th);
else
    bubble_out = repetitive_ellipse(bubble,maxsize_th,minsize_th,ratio_th,maxratio_th);
end

end


function bubble_out = fitting_del(bubble,fitting_th)
[num2,~] = size(bubble);
%删除与拟合边界拟合度过小的椭圆
for n=1:num2
    if(abs(bubble{n,2})>1)
        [b,~]=size(bubble{n,6});
        for i=1:b
            if(bubble{n,6}{i,1}~=0)
                fitting_rate = my_arc_fitting(bubble{n,7}{i,1},'method_ellipse',...
                            [bubble{n,6}{i,1}(1,3),bubble{n,6}{i,1}(1,4),bubble{n,6}{i,1}(1,1),bubble{n,6}{i,1}(1,2),bubble{n,6}{i,1}(1,5)]);
                if(fitting_rate < fitting_th)
                    bubble{n,6}{i,1}=0;
                    warning(['气泡',num2str(n),'第',num2str(i),'个在拟合度不过关被删除，其拟合度为：',num2str(fitting_rate)]);
                end
            end
        end
    end
    if bubble{n,2} == 0 
        fitting_rate = my_arc_fitting(bubble{n,7},'method_ellipse',...
                    [bubble{n,6}(1,3),bubble{n,6}(1,4),bubble{n,6}(1,1),bubble{n,6}(1,2),bubble{n,6}(1,5)]);
        if(fitting_rate < fitting_th)
            bubble{n,6}=0;
            warning(['气泡',num2str(n),'在拟合度不过关被删除，其拟合度为：',num2str(fitting_rate)]);
        end
    end
end

bubble_out = bubble;
end


function bubble_out = abnormal_ellipse(bubble,maxsize_th,minsize_th,ratio_th,maxratio_th)
[num2,~] = size(bubble);

% 剔除异常椭圆
for n=1:num2
    if(abs(bubble{n,2})>1)
        [b,~]=size(bubble{n,6});
        for i=1:b
            if(bubble{n,6}{i,1}~=0)
                if(bubble{n,6}{i,1}(1,1)<minsize_th||bubble{n,6}{i,1}(1,2)<minsize_th||bubble{n,6}{i,1}(1,1)>maxsize_th|| ...
                        bubble{n,6}{i,1}(1,2)>maxsize_th||bubble{n,6}{i,1}(1,1)/bubble{n,6}{i,1}(1,2)<0.4||bubble{n,6}{i,1}(1,1)/bubble{n,6}{i,1}(1,2)>maxratio_th)
                    bubble{n,6}{i,1}=0;
                    warning(['气泡',num2str(n),'第',num2str(i),'个在气泡异常中被删除']);
                end
            end
        end
    end
    if bubble{n,2} == 1 & bubble{n,6}~=0
        if(bubble{n,6}(1,1)<minsize_th||bubble{n,6}(1,2)<minsize_th||bubble{n,6}(1,1)>maxsize_th|| ...
                        bubble{n,6}(1,2)>maxsize_th || bubble{n,6}(1,1)/bubble{n,6}(1,2)<ratio_th || bubble{n,6}(1,1)/bubble{n,6}(1,2)>maxratio_th)
            bubble{n,6}=0;
            warning(['气泡',num2str(n),'在气泡异常中被删除']);
        end
    end
    if bubble{n,2} == 0 & bubble{n,6}~=0
        if(bubble{n,6}(1,1)>maxsize_th || bubble{n,6}(1,2)>maxsize_th || bubble{n,6}(1,1)/bubble{n,6}(1,2)<ratio_th || bubble{n,6}(1,1)/bubble{n,6}(1,2)>maxratio_th)
            bubble{n,6}=0;
            warning(['气泡',num2str(n),'在气泡异常中被删除']);
        end
    end
end

bubble_out = bubble;
end

function bubble_out = repetitive_ellipse(bubble,repetition_length_th)
[num2,~] = size(bubble);
%剔除重复椭圆
for n=1:num2
    if(abs(bubble{n,2})>1)
        [b,~]=size(bubble{n,6});
        for i=1:(b-1)
            for j=(i+1):b
                if(isempty(bubble{n,6}{i,1})~=1 & isempty(bubble{n,6}{j,1})~=1 & bubble{n,6}{i,1}~=0 & bubble{n,6}{j,1}~=0)
                    if(sqrt((bubble{n,6}{i,1}(1,3)-bubble{n,6}{j,1}(1,3))^2+(bubble{n,6}{i,1}(1,4)-bubble{n,6}{j,1}(1,4))^2) < repetition_length_th)
                        % 以下为2021、6、28 添加 增加重复椭圆的判定准则
                        fitting_rate_1 = my_arc_fitting(bubble{n,7}{i,1},'method_ellipse',...
                            [bubble{n,6}{i,1}(1,3),bubble{n,6}{i,1}(1,4),bubble{n,6}{i,1}(1,1),bubble{n,6}{i,1}(1,2),bubble{n,6}{i,1}(1,5)]);
                        fitting_rate_2 = my_arc_fitting(bubble{n,7}{j,1},'method_ellipse',...
                            [bubble{n,6}{j,1}(1,3),bubble{n,6}{j,1}(1,4),bubble{n,6}{j,1}(1,1),bubble{n,6}{j,1}(1,2),bubble{n,6}{j,1}(1,5)]);
                        
                        if fitting_rate_1 > fitting_rate_2
                            bubble{n,6}{j,1}=0;
                            warning(['气泡',num2str(n),'第',num2str(j),'个在气泡重复中被删除，于其重复的是气泡',num2str(n),'第',num2str(i),'拟合度为',num2str(fitting_rate_1)]);
                        else
                            bubble{n,6}{i,1}=0;
                            warning(['气泡',num2str(n),'第',num2str(i),'个在气泡重复中被删除，于其重复的是气泡',num2str(n),'第',num2str(j),'拟合度为',num2str(fitting_rate_2)]);
                        end
                        % 原来 bubble{n,6}{j,1}=0;
                    end
                end
            end
        end
    end
end

bubble_out = bubble;
end
