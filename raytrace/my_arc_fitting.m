function fitting_rate = my_arc_fitting(arc,varargin)
% 本函数的作用是为气泡图像处理算法提供分段弧与参考边界的拟合度计算，使用方式如下：
% 
% 输入：
%     arc：弧
%     reference_bw：参考边界
%     varargin：可变长度输入列表，用于指定阈值
% 
% 输出：
%     fitting_rate：拟合度计算结果
% 
% 调用说明：
%     fitting_rate = my_arc_fitting(arc,reference_bw)
%                  默认预处理方式，输出为拟合度计算结果
%     fitting_rate = my_arc_fitting(arc,reference_bw,5)
%                  设定弧点拟合阈值为5像素，即弧点与参考边界的最小距离阈值为5
% 
% 版本号VOL1.0，编写于2021年6月28日，作者：WG-Chen
% 版本号VOL2.0，编写于2021年6月29日，作者：WG-Chen
% 修改内容：增加参考边界为离散椭圆边界的方法，具体方案是先做出椭圆，找到椭圆上离弧段一点距离最小的两个点；
%          然后在两个点中加密弧段长度个点，再寻找弧段一点离椭圆的最小距离；


%% default value 
l_th = 3;
method = 1;
arc_num = 1e2;

%% methods

if numel(varargin) == 1 
    reference_bw = varargin{1};
elseif numel(varargin) > 1
    tf = {'length_th','method_ellipse','arc_num'};
    for i = 1 : numel(tf)
        for j = 1 : numel(varargin)
            tf_1 = strcmp(tf{i},varargin{j});
            if tf_1 == 1 
                switch i
                    case 1
                        l_th = varargin{j+1};
                    case 2
                        method = 2;
                        ell_data = varargin{j+1};
                        TF = isnan(ell_data);
                        TF_sum = sum(TF);
                    case 3
                        arc_num = varargin{j+1};
                end
            end
        end
    end
end


%% 计算部

if method == 1   % 参考边界法
    count = 0;

    [arc_l,~] = size(arc);

    for i = 1 : arc_l
        leng = zeros();
        [re_l,~] = size(reference_bw);
        for j = 1 : re_l
            leng(j) = norm(arc(i,:)-reference_bw(j,:));
        end
        min_len = min(leng);
        A = leng == min_len;
        if min_len <= l_th
            count = count + 1 ;
            reference_bw(A,:) = [];
        end
        clear length re_l j min_len A
    end

    fitting_rate = count/arc_l;
    
else % 椭圆边界法
    if TF_sum > 0 
        fitting_rate = 0;
    else
        count = 0;
        [arc_l,~] = size(arc);
        if arc_l < 100
            arc_num = arc_l;
        end
        elli_1 = write_ellipse(ell_data(1),ell_data(2),ell_data(3),ell_data(4),ell_data(5),arc_num);

        % first try
        for i = 1 : arc_l
            leng = zeros();
            [re_l,~] = size(elli_1);
            for j = 1 : re_l
                leng(j) = norm(arc(i,:)-elli_1(j,:));
            end
            min_len = min(leng);
            A_loc = find(leng == min_len);
            A_loc = A_loc(1);
            leng(A_loc) = 1000000;
            min_len = min(leng);
            B_loc = find(leng == min_len);
            B_loc = B_loc(1);


            % second try
            t_2 = linspace(1,arc_num+1,arc_num);
            xita_1 = 2*pi/arc_num*t_2(A_loc);
            xita_2 = 2*pi/arc_num*t_2(B_loc);
            xita = abs(xita_1 - xita_2);
            if xita_1 <= xita_2
                xita_begin = xita_1;
            else
                xita_begin = xita_2;
            end
            x_2 = ell_data(1)+ell_data(3)*cos(xita_begin + xita/arc_num*t_2)*cos(ell_data(5))-ell_data(4)*sin(xita_begin + xita/arc_num*t_2)*sin(ell_data(5));
            y_2 = ell_data(2)-ell_data(3)*cos(xita_begin + xita/arc_num*t_2)*sin(ell_data(5))-ell_data(4)*sin(xita_begin + xita/arc_num*t_2)*cos(ell_data(5));
            elli_2 = zeros(length(x_2),2);
            elli_2(:,1) = y_2;
            elli_2(:,2) = x_2;   % 2021/6/29 修改，注意这里x，y的顺序！！！

            leng = zeros();
            [re_l,~] = size(elli_2);
            for j = 1 : re_l
                leng(j) = norm(arc(i,:)-elli_2(j,:));
            end
            min_len = min(leng);
            if min_len <= l_th
                count = count + 1 ;
            end

        end

        fitting_rate = count/arc_l;
    end

end

end


% 2021/6/28 添加
function elli = write_ellipse(x0,y0,a,b,theta,num_t)
%%  说明
%本程序画一个中心在（x0，y0)处的椭圆，其长短轴分别为a,b,椭圆沿Z轴转theta角  的所有点
%%
t = linspace(1,num_t+1,num_t);
x = x0+a*cos(2*pi/num_t*t)*cos(theta)-b*sin(2*pi/num_t*t)*sin(theta);
y = y0-a*cos(2*pi/num_t*t)*sin(theta)-b*sin(2*pi/num_t*t)*cos(theta);
elli = zeros(length(x),2);
elli(:,1) = y;
elli(:,2) = x;
end

