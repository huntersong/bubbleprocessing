function [bubble,bub_overlap,ao_data] = my_bub_processing(bub_bw,lab_image,varargin)
% 本函数的作用是为气泡图像处理算法提供气泡边界处理功能，对单气泡进行分类拟合处理，使用方式如下：
% 
% 输入：
%     bub_bw：气泡边界信息
%     lab_image：带标记的气泡区域
%     varargin：可变长度输入列表，用于指定阈值
% 
% 输出：
%     bubble：识别到的气泡信息，包括圆度、质心、尺寸、角度等，用于气泡后处理my_overlapbubbles
%     bub_overlap：重叠气泡边界数组，用于调用重叠气泡处理模块my_overlapbubbles
%     ao_data：重叠气泡的凹点信息
% 
% 调用说明：
%     bubble = my_bub_processing(bub_bw,lab_image,varargin)
%                  默认预处理方式，输出为识别的气泡信息
%     bubble = my_bub_processing(bub_bw,lab_image,'roundness',0.94,'point-apart',10,'bubble-size-th',40);
%                  可识别关键词：  'roundness'（圆度阈值） 'point-apart'（间隔阈值）
%                  'bubble-size-th' （凹点识别边界尺寸阈值）
%                  'length_width_ratio'（长宽比阈值） 'method'(0 圆度 1 长宽比 2 综合判定)
%                  'single_ao_length'（单凹点边界去除阈值） 
% 
% 版本号VOL1.0，编写于2021年6月3日，作者：WG-Chen

%% default value  
roundness_th = 0.94;                % 圆度默认判定阈值
point_apart = 5;                    % 前序点以及后续点的默认间隔 unchangable!
bub_size_th = 30;                   % 凹点识别的默认尺寸阈值，   ！不允许小于2倍的单凹点判定阈值！
length_width_ratio = [0.95,1.1];    % 长宽比默认判定阈值   unchangable
method_flag = 1;                    % 判定圆、椭圆拟合的方法；默认值为1：长宽比判定法；圆度法：0；同时采用两种方法：2； unchangbale
single_ao_length = 10;              % 单凹点去除两侧的边界阈值 unchangable
bubble_boundaries = bub_bw;
[bub_num,~] = size(bubble_boundaries);

%% methods
if numel(varargin) > 0
    tf = {'roundness','point-apart','bubble-size-th','length_width_ratio','method','single_ao_length'};
    for i = 1 : numel(tf)
        for j = 1 : numel(varargin)
            tf = strcmp(tf{i},varargin{j});
            if tf == 1
                switch i
                    case 1
                        roundness_th = varargin{j+1};
                    case 2
                        point_apart = varargin{j+1};
                    case 3
                        bub_size_th = varargin{j+1};
                    case 4
                        length_width_ratio = varargin{j+1};
                    case 5
                        method_flag = varargin{j+1};
                    case 6
                        single_ao_length = varargin{j+1};  
                end
            end
        end
    end
end
%% 计算部
bubble_boundaries_l = zeros();
for i = 1 : bub_num
    [bubble_boundaries_l(i,1),~] = size(bubble_boundaries{i,1});
end
bubble_boundary_number=max(bubble_boundaries_l(:,1))+1;                     %气泡边界识别列数

ii=zeros(bub_num,bubble_boundary_number);                                      %气泡边界行数
jj=zeros(bub_num,bubble_boundary_number);                                      %气泡边界列数
 for n=1:bub_num                                                               %从上到下，从左往右，识别气泡第一个点
     lx = length(bubble_boundaries{n}(:,1));
     ii(n,1:lx-1) = bubble_boundaries{n,1}(1:lx-1,1);
     jj(n,1:lx-1) = bubble_boundaries{n,1}(1:lx-1,2);
 end

iii=zeros(bub_num,bubble_boundary_number);                                     %过滤气泡边界过短的气泡
jjj=zeros(bub_num,bubble_boundary_number);
for n=1:bub_num
    if(ii(n,2)==0||ii(n,1)==0)
        continue
    else
        iii(n,:)=ii(n,:);
        jjj(n,:)=jj(n,:);
    end
end

F_bubble_num=zeros(bub_num,1);                                                 %气泡边界长度
for n=1:bub_num 
    for k=1:bubble_boundary_number
        if(iii(n,k)~=0&&iii(n,k+1)==0)
            F_bubble_num(n,1)=k;
        end
    end
end
X=zeros(bub_num,2);                                                            %气泡周围绿框的xy值
Y=zeros(bub_num,2);
for n=1:bub_num
    if(iii(n,1)>0)
        X(n,1)=min(jjj(n,1:F_bubble_num(n,1)));
        X(n,2)=max(jjj(n,1:F_bubble_num(n,1)));
        Y(n,1)=min(iii(n,1:F_bubble_num(n,1)));
        Y(n,2)=max(iii(n,1:F_bubble_num(n,1)));
    end
end

%% 识别凹点
%针对气泡边界数量大于30的气泡进行气泡识别
%三角函数识别
theta=zeros(bub_num,bubble_boundary_number);                                 %通过坐标计算转动角度的方法（不够准确）
aa1=zeros(bub_num,bubble_boundary_number);                                   %后续点角度
aa2=zeros(bub_num,bubble_boundary_number);                                   %前序点角度
aaa=zeros(bub_num,bubble_boundary_number);                                   %夹角
for n=1:bub_num
    for k=1+point_apart:bubble_boundary_number-point_apart
        if(iii(n,k)==0)
            break
        elseif(iii(n,k+point_apart)==0)
            break
        elseif(F_bubble_num(n,1) < bub_size_th)                           %气泡周长小于30像素的忽略
            break
        else
            x1=jjj(n,k-point_apart);                                      %后续点
            y1=iii(n,k-point_apart);                                      %后续点
            x2=jjj(n,k);
            y2=iii(n,k);
            x3=jjj(n,k+point_apart);                                      %前序点
            y3=iii(n,k+point_apart);                                      %前序点
            theta(n,k)=acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x1-x2,y3-y2])));
            if x2~=x1
                k1=(y2-y1)/(x1-x2);                                       %斜率
            end
            % 夹角等于前序点逆时针到x轴正半轴的角度加上后续点顺时针到x轴正半轴的角度和
            if x3~=x2
                k2=(y3-y2)/(x2-x3);
            end
            %后续点采用从x轴正半轴开始逆时针所需要转动角度
            if(x1>x2)
                if(y1>y2) %第四象限
                    a1=360-atand(abs(k1));
                elseif(y1<y2) %第一象限
                    a1=atand(k1);
                elseif(y1==y2) %x正轴
                    a1=0;
                end
            elseif(x1<x2) %第二三象限
                a1=180+atand(k1);
            elseif(x1==x2)
                if(y1>y2)
                    a1=270;%y轴负半轴
                else
                    a1=90;%y轴正半轴
                end
            end

            if(x3>x2)
                if(y3<y2)   %第一象限
                    a2=atand(k2);
                elseif(y3>y2)  %第四象限
                    a2=360-atand(abs(k2));
                elseif(y3==y2) %X轴正半轴
                    a2=0;
                end
            elseif(x3<x2)
                a2=180+atand(k2);
            elseif(x3==x2)
                if(y3>y2)
                    a2=270;  %y轴负半轴
                else
                    a2=90;   %y轴正半轴
                end
            end
            aa1(n,k)=a1;
            aa2(n,k)=a2;
        end
    end
end
%计算前序点逆时针转动到后续点的角度，分类讨论
for n=1:bub_num
    for k=1+point_apart:bubble_boundary_number-point_apart
        if(iii(n,k)==0)
            break
        elseif(iii(n,k+point_apart)==0)
            break
        else
            if(aa2(n,k)>aa1(n,k))
                aaa(n,k)=aa2(n,k)-aa1(n,k);
            elseif(aa1(n,k)>aa2(n,k))
                aaa(n,k)=360-(aa1(n,k)-aa2(n,k));
            end
        end
    end
end

aaaa=zeros(bub_num,bubble_boundary_number);                                  %角度过滤
for n=1:bub_num
    if(F_bubble_num(n,1)>11)
        aaaa(n,point_apart+1)=0.8*aaa(n,point_apart+1)+0.2*aaa(n,point_apart+2);
        aaaa(n,F_bubble_num(n,1)-point_apart)=0.8*aaa(n,F_bubble_num(n,1)-point_apart)+0.2*aaa(n,F_bubble_num(n,1)-point_apart-1);
        for i=point_apart+2:F_bubble_num(n,1)-point_apart-1
            aaaa(n,i)=0.2*aaa(n,i-1)+0.6*aaa(n,i)+0.2*aaa(n,i+1);
        end
    end
end

%%  识别极大值点
%相邻两点求梯度，梯度最大区域进行匹配。第一个上升值 up ，第一个下降值 down 进行匹配。在up(1,1) down(1,1)区间内选取最大值。
%首先判断是否存在大于210度的点，在大于210度的点处进行极大值分析
ao_pks=zeros(bub_num,bubble_boundary_number);
ao_locs=zeros(bub_num,bubble_boundary_number);
for n=1:bub_num     %获取最大
    if(max(aaaa(n,:))<180)
        continue
    else
        [pks,locs] = findpeaks(aaaa(n,:));
        ao_pks(n,2:length(pks)+1)=pks;
        ao_locs(n,2:length(locs)+1)=locs;
    end
end

for n=1:bub_num
    counter = ao_pks(n,:)>1;
    for i=1:sum(double(counter))                                            %判断峰值周围是否有更大的峰值
        if (ao_locs(n,i)>10&&ao_locs(n,i)<(F_bubble_num(n,1)-10))
            if(sum(aaaa(n,ao_locs(n,i))<aaaa(n,ao_locs(n,i)-10:ao_locs(n,i)+10))>0)
                ao_locs(n,i)=0;
                ao_pks(n,i)=0;
            end
        end
    end
end

ao_point_num=zeros(bub_num,1);                                               %凹点数量
ao_point_de=zeros(bub_num,bubble_boundary_number);                           %凹点的在ao_locs中的位置 ao_locs(n,(ao_point_de(n,i)))
ao_point_de2=cell(bub_num,bubble_boundary_number);                           %凹点的坐标值  jjj(n,ao_locs(n,(ao_point_de(n,i))))
theta_th=200;                                                             %凹点判断阈值
for n=1:bub_num
    ao_point_num(n,1)=length(find(ao_pks(n,:)>theta_th));
    ao_point_de(n,1:length(find(ao_pks(n,:)>theta_th)))=find(ao_pks(n,:)>theta_th);
    if(ao_point_num(n,1)~=0)
        for i=1:ao_point_num(n,1)
            ao_point_de2{n,i}=[int16(iii(n,ao_locs(n,(ao_point_de(n,i))))),int16(jjj(n,ao_locs(n,(ao_point_de(n,i)))))];%凹点坐标
        end
    end
end

%% 识别信息输出
% 1.是否斑点 2.凹点个数 3.圆度 4.面积 5.长宽比 6.气泡拟合信息
bubble = cell(bub_num,6);
for n=1:bub_num
    if(iii(n,1)==0)
        bubble{n,1} = 1;
    else
        bubble{n,2} = ao_point_num(n,1);
        bubble{n,1} = 0;
    end
end
for n=1:bub_num
    s = length(find(lab_image==bubble_boundaries{n,2}));
    bubble{n,4}= s;
    bubble{n,3}= 4*pi*bubble{n,4}/(F_bubble_num(n,1)^2);
    bubble{n,5}= (X(n,2)-X(n,1))/(Y(n,1)-Y(n,2));
end

%根据长宽比 决定圆拟合还是椭圆拟合
% 无凹点类型的
if method_flag == 1
    for n=1:bub_num
        if(abs(bubble{n,2})==0)
            if(abs(bubble{n,5})>length_width_ratio(1) && abs(bubble{n,5})<length_width_ratio(2))  %长宽比判定准则
                [radius,cir_xc,cir_yc] = my_circfit(iii(n,1:F_bubble_num(n,1)),jjj(n,1:F_bubble_num(n,1)));
                bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,6}(1,1) = radius;
                bubble{n,6}(2,1) = round(cir_xc);
                bubble{n,6}(2,2) = round(cir_yc);
                bubble{n,2}=0.5;    % 无凹点圆拟合
            else
                [ bubble{n,6}(1,1), bubble{n,6}(1,2), bubble{n,6}(1,3), ...
                    bubble{n,6}(1,4), bubble{n,6}(1,5)] =  ...
                    my_ellipsefit(jjj(n,1:F_bubble_num(n,1)),iii(n,1:F_bubble_num(n,1)));
                    bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                    bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,2}=0;      % 无凹点椭圆拟合
            end
        end
    end
elseif method_flag == 0
    for n=1:bub_num
        if(abs(bubble{n,2})==0)
            if(abs(bubble{n,3})>roundness_th)  %圆度判定准则
                [radius,cir_xc,cir_yc] = my_circfit(iii(n,1:F_bubble_num(n,1)),jjj(n,1:F_bubble_num(n,1)));
                bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,6}(1,1) = radius;
                bubble{n,6}(2,1) = round(cir_xc);
                bubble{n,6}(2,2) = round(cir_yc);
                bubble{n,2}=0.5;    % 无凹点圆拟合
            else
                [ bubble{n,6}(1,1), bubble{n,6}(1,2), bubble{n,6}(1,3), ...
                    bubble{n,6}(1,4), bubble{n,6}(1,5)] =  ...
                    my_ellipsefit(jjj(n,1:F_bubble_num(n,1)),iii(n,1:F_bubble_num(n,1)));
                    bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                    bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,2}=0;      % 无凹点椭圆拟合
            end
        end
    end
else
    for n=1:bub_num
        if(abs(bubble{n,2})==0)
            if(abs(bubble{n,3})> roundness_th && abs(bubble{n,5}) > length_width_ratio(1) && abs(bubble{n,5}) < length_width_ratio(2))  %综合判定准则
                [radius,cir_xc,cir_yc] = my_circfit(iii(n,1:F_bubble_num(n,1)),jjj(n,1:F_bubble_num(n,1)));
                bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,6}(1,1) = radius;
                bubble{n,6}(2,1) = round(cir_xc);
                bubble{n,6}(2,2) = round(cir_yc);
                bubble{n,2}=0.5;    % 无凹点圆拟合
            else
                [ bubble{n,6}(1,1), bubble{n,6}(1,2), bubble{n,6}(1,3), ...
                    bubble{n,6}(1,4), bubble{n,6}(1,5)] =  ...
                    my_ellipsefit(jjj(n,1:F_bubble_num(n,1)),iii(n,1:F_bubble_num(n,1)));
                    bubble{n,7}(:,1) = iii(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                    bubble{n,7}(:,2) = jjj(n,1:F_bubble_num(n,1));    % 2021/06/28 添加
                bubble{n,2}=0;      % 无凹点椭圆拟合
            end
        end
    end
end

%单凹点 边界去除凹点两侧20个单位后进行椭圆拟合
single_ao_i=zeros(bub_num,bubble_boundary_number);
single_ao_j=zeros(bub_num,bubble_boundary_number);

for n=1:bub_num
    if(abs(bubble{n,2})==1)
        if(ao_locs(n,(ao_point_de(n,1)))>single_ao_length&&(ao_locs(n,(ao_point_de(n,1)))+single_ao_length)<F_bubble_num(n,1))
            single_ao_i(n,1:F_bubble_num(n,1)-single_ao_length-ao_locs(n,(ao_point_de(n,1))))=iii(n,ao_locs(n,(ao_point_de(n,1)))+single_ao_length+1:F_bubble_num(n,1));
            single_ao_i(n,F_bubble_num(n,1)-single_ao_length-ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-2*single_ao_length)=iii(n,1:ao_locs(n,(ao_point_de(n,1)))-single_ao_length);
            single_ao_j(n,1:F_bubble_num(n,1)-single_ao_length-ao_locs(n,(ao_point_de(n,1))))=jjj(n,ao_locs(n,(ao_point_de(n,1)))+single_ao_length+1:F_bubble_num(n,1));
            single_ao_j(n,F_bubble_num(n,1)-single_ao_length-ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-2*single_ao_length)=jjj(n,1:ao_locs(n,(ao_point_de(n,1)))-single_ao_length);
        elseif(ao_locs(n,(ao_point_de(n,1)))<=single_ao_length)
            single_ao_i(n,1:F_bubble_num(n,1)-2*single_ao_length)=iii(n,ao_locs(n,(ao_point_de(n,1)))+single_ao_length+1:F_bubble_num(n,1)-(single_ao_length-ao_locs(n,(ao_point_de(n,1)))));
            single_ao_j(n,1:F_bubble_num(n,1)-2*single_ao_length)=jjj(n,ao_locs(n,(ao_point_de(n,1)))+single_ao_length+1:F_bubble_num(n,1)-(single_ao_length-ao_locs(n,(ao_point_de(n,1)))));
        elseif(ao_locs(n,(ao_point_de(n,1)))>=F_bubble_num(n,1)-single_ao_length)
            single_ao_i(n,1:F_bubble_num(n,1)-2*single_ao_length)=iii(n,(single_ao_length-(F_bubble_num(n,1)-ao_locs(n,(ao_point_de(n,1)))))+1:ao_locs(n,(ao_point_de(n,1)))-single_ao_length);
            single_ao_j(n,1:F_bubble_num(n,1)-2*single_ao_length)=jjj(n,(single_ao_length-(F_bubble_num(n,1)-ao_locs(n,(ao_point_de(n,1)))))+1:ao_locs(n,(ao_point_de(n,1)))-single_ao_length);
        end
    end
end

for n=1:bub_num
    if(abs(bubble{n,2})==1)
        [ bubble{n,6}(1,1), bubble{n,6}(1,2), bubble{n,6}(1,3), bubble{n,6}(1,4), bubble{n,6}(1,5)] =  ...
            my_ellipsefit(single_ao_j(n,1:F_bubble_num(n,1)-2*single_ao_length),single_ao_i(n,1:F_bubble_num(n,1)-2*single_ao_length));
        bubble{n,7}(:,1) = single_ao_i(n,1:F_bubble_num(n,1)-2*single_ao_length);    % 2021/06/28 添加
        bubble{n,7}(:,2) = single_ao_j(n,1:F_bubble_num(n,1)-2*single_ao_length);    % 2021/06/28 添加
    end
end

%双凹点
%将气泡范围分区
double_ao_i=cell(bub_num,2);
double_ao_j=cell(bub_num,2);
for n=1:bub_num
    if(abs(bubble{n,2})==2)
        double_ao_i{n,1} = iii(n,ao_locs(n,(ao_point_de(n,1)))+1:ao_locs(n,(ao_point_de(n,2))));%两凹点顺时针边界
        double_ao_i{n,2}(1,1:ao_locs(n,(ao_point_de(n,1)))) = iii(n,ao_locs(n,(ao_point_de(n,1))):-1:1);%第二段逆时针边界
        double_ao_i{n,2}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,2))))+ao_locs(n,(ao_point_de(n,1)))) = ...
            iii(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,2)))+1);
        double_ao_j{n,1} = jjj(n,ao_locs(n,(ao_point_de(n,1)))+1:ao_locs(n,(ao_point_de(n,2))));
        double_ao_j{n,2}(1,1:ao_locs(n,(ao_point_de(n,1)))) = jjj(n,ao_locs(n,(ao_point_de(n,1))):-1:1);%第二段逆时针边界
        double_ao_j{n,2}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,2))))+ao_locs(n,(ao_point_de(n,1)))) = ...
            jjj(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,2)))+1);
    end
end
for n=1:bub_num
    if(abs(bubble{n,2})==2)
        [ bubble{n,6}{1,1}(1,1), bubble{n,6}{1,1}(1,2), bubble{n,6}{1,1}(1,3),bubble{n,6}{1,1}(1,4), bubble{n,6}{1,1}(1,5)] =  ...
            my_ellipsefit(double_ao_j{n,1}(1,1:length(double_ao_j{n,1})),double_ao_i{n,1}(1,1:length(double_ao_i{n,1})));
        bubble{n,7}{1,1}(:,1) = double_ao_i{n,1}(1,1:length(double_ao_i{n,1}));    % 2021/06/28 添加
        bubble{n,7}{1,1}(:,2) = double_ao_j{n,1}(1,1:length(double_ao_j{n,1}));    % 2021/06/28 添加
        [ bubble{n,6}{2,1}(1,1), bubble{n,6}{2,1}(1,2), bubble{n,6}{2,1}(1,3),bubble{n,6}{2,1}(1,4), bubble{n,6}{2,1}(1,5)] =  ...
            my_ellipsefit(double_ao_j{n,2}(1,1:length(double_ao_j{n,2})),double_ao_i{n,2}(1,1:length(double_ao_i{n,2})));
        bubble{n,7}{2,1}(:,1) = double_ao_i{n,2}(1,1:length(double_ao_i{n,2}));    % 2021/06/28 添加
        bubble{n,7}{2,1}(:,2) = double_ao_j{n,2}(1,1:length(double_ao_j{n,2}));    % 2021/06/28 添加
    end
end


% 至此对给定的气泡边界进行了完全分类，并对无凹点进行判定后拟合，对单凹点进行删除凹点附近边界后拟合，
% 对双凹点进行分区后拟合，输出为已识别、分类的气泡边界、重叠气泡信息
% output:

bub_overlap = {};
overlap_num = 0;
for n = 1 : bub_num
    if bubble{n,2} > 2
        overlap_num = overlap_num + 1;
        bub_overlap{overlap_num,1} = bubble{n,1};
        bub_overlap{overlap_num,2} = bubble{n,2};
        bub_overlap{overlap_num,3} = bubble{n,3};
        bub_overlap{overlap_num,4} = bubble{n,4};
        bub_overlap{overlap_num,5} = bubble{n,5};
        bub_overlap{overlap_num,6} = bubble{n,6};
        bub_overlap{overlap_num,7} = n;              % 记录重叠气泡在原气泡边界数组的位置
        mult_ao_i = {};mult_ao_j = {};
        for i=1:abs(bubble{n,2})-1
            mult_ao_i{1,i}=iii(n,ao_locs(n,(ao_point_de(n,i)))+1:ao_locs(n,(ao_point_de(n,i+1))));
            mult_ao_j{1,i}=jjj(n,ao_locs(n,(ao_point_de(n,i)))+1:ao_locs(n,(ao_point_de(n,i+1))));
        end
        mult_ao_i{1,abs(bubble{n,2})}(1,1:ao_locs(n,(ao_point_de(n,1))))=iii(n,ao_locs(n,(ao_point_de(n,1))):-1:1);
        mult_ao_i{1,abs(bubble{n,2})}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,abs(bubble{n,2})))))+ao_locs(n,(ao_point_de(n,1))))...
            =iii(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,abs(bubble{n,2}))))+1);
        mult_ao_j{1,abs(bubble{n,2})}(1,1:ao_locs(n,(ao_point_de(n,1))))=jjj(n,ao_locs(n,(ao_point_de(n,1))):-1:1);
        mult_ao_j{1,abs(bubble{n,2})}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,abs(bubble{n,2})))))+ao_locs(n,(ao_point_de(n,1))))...
            =jjj(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,abs(bubble{n,2}))))+1);
        mult_ao_ij = cell(1,1);
        for ha = 1 : length(mult_ao_i)
            mult_ao_ij{1,1}{ha,1}(:,1) = mult_ao_i{1,ha}(1,:);
            mult_ao_ij{1,1}{ha,1}(:,2) = mult_ao_j{1,ha}(1,:);
        end
        bub_overlap{overlap_num,8} = mult_ao_ij;
    end
end


ao_data = cell(overlap_num,6); % ao_locs,ao_point_de,ao_point_de2,iii,jjj,F_bubble_num
for n = 1 : overlap_num
    ao_data{n,1} = ao_locs(bub_overlap{n,7},:);
    ao_data{n,2} = ao_point_de(bub_overlap{n,7},:);
    for nn = 1 : bubble_boundary_number
        ao_data{n,3}{1,nn} = ao_point_de2{bub_overlap{n,7},nn};
    end
    ao_data{n,4} = iii(bub_overlap{n,7},:);
    ao_data{n,5} = jjj(bub_overlap{n,7},:);
    ao_data{n,6} = F_bubble_num(bub_overlap{n,7},:);
end


end