function bub_overlap = my_overlapbubbles(bubble,ao_data,image,varargin)
% 本函数的作用是为气泡图像处理算法提供气泡边界处理功能，对单气泡进行分类拟合处理，使用方式如下：
% 
% 输入：
%     bubble：重叠气泡边界信息
%     N×6的元胞数组，{1}存储是否为斑点气泡；{2}存储凹点个数；{3}为圆度；{4}为面积；{5}为长宽比；{6}储存为拟合信息
%     ao_data：重叠气泡凹点信息
%     image：去背景后的原始图像
%     varargin：可变长度输入列表，用于指定阈值
% 
% 输出：
%     bub_overlap：重叠气泡识别结果
% 
% 调用说明：
%     bub_overlap = my_overlapbubbles(bubble)
%                  默认预处理方式，输出为识别的重叠气泡信息
%     bub_overlap = my_overlapbubbles(bubble,'major_arc_th',0.94,);
%                  设定圆度为0.94
% 
% 版本号VOL1.0，编写于2021年6月3日，作者：WG-Chen
%% default value
major_arc_th = 25; %优弧阈值
partner_cluster_distance = 0.05;
partner_cluster_proportion = 0.95;
major_arc_angle = 85;

%% methods
if nargin > 3
    major_arc_th = varargin{2};
end

%% input:
[num2,~] = size(bubble); 
[~,bubble_boundary_number] = size(ao_data{1,3});
for n = 1 : num2
    ao_locs(n,:) = ao_data{n,1};
    ao_point_de(n,:) = ao_data{n,2};
    for nn = 1 : bubble_boundary_number
        ao_point_de2{n,nn} = ao_data{n,3}{1,nn};
    end
    iii(n,:) = ao_data{n,4};
    jjj(n,:) = ao_data{n,5};
    F_bubble_num(n,:) = ao_data{n,6};
    ao_point_num(n,1) = bubble{n,2};
end

%% 计算部

%多凹点
mult=max(ao_point_num);
mult_ao_i=cell(num2,mult);
mult_ao_j=cell(num2,mult);
mult_finish=ones(num2,1);

%多凹点分弧
for n=1:num2
    for i=1:abs(bubble{n,2})-1
        mult_ao_i{n,i}=iii(n,ao_locs(n,(ao_point_de(n,i)))+1:ao_locs(n,(ao_point_de(n,i+1))));
        mult_ao_j{n,i}=jjj(n,ao_locs(n,(ao_point_de(n,i)))+1:ao_locs(n,(ao_point_de(n,i+1))));
    end
    mult_ao_i{n,abs(bubble{n,2})}(1,1:ao_locs(n,(ao_point_de(n,1))))=iii(n,ao_locs(n,(ao_point_de(n,1))):-1:1);
    mult_ao_i{n,abs(bubble{n,2})}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,abs(bubble{n,2})))))+ao_locs(n,(ao_point_de(n,1))))...
        =iii(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,abs(bubble{n,2}))))+1);
    mult_ao_j{n,abs(bubble{n,2})}(1,1:ao_locs(n,(ao_point_de(n,1))))=jjj(n,ao_locs(n,(ao_point_de(n,1))):-1:1);
    mult_ao_j{n,abs(bubble{n,2})}(1,ao_locs(n,(ao_point_de(n,1)))+1:F_bubble_num(n,1)-(ao_locs(n,(ao_point_de(n,abs(bubble{n,2})))))+ao_locs(n,(ao_point_de(n,1))))...
        =jjj(n,F_bubble_num(n,1):-1:ao_locs(n,(ao_point_de(n,abs(bubble{n,2}))))+1);
end

%每段弧中点
mult_mid_point=cell(num2,mult);
for n=1:num2
    for i=1:abs(bubble{n,2})
        mult_mid_point{n,i}(1,1)=round(length(mult_ao_i{n,i})/2);
        mult_mid_point{n,i}(2,1)=mult_ao_i{n,i}(1,mult_mid_point{n,i}(1,1));
        mult_mid_point{n,i}(2,2)=mult_ao_j{n,i}(1,mult_mid_point{n,i}(1,1));
    end
end

%%
%判断每段弧的情况
mult_flag=cell(num2,mult);      %多凹点气泡每个弧对应的情况 (1,1)是否是优弧(1,2)类型(2,1)是否拟合(3,1)对应弧段
dis_A_mid=zeros(num2,mult);
dis_B_mid=zeros(num2,mult);
dis_A_B=zeros(num2,mult);
%记录所有弧段的横纵坐标
extreme_point=cell(num2,20);
for n=1:num2
    for i=1:abs(bubble{n,2})
        extreme_point{n,i}(1,1)=min(mult_ao_i{n,i}); %i_min
        extreme_point{n,i}(1,2)=max(mult_ao_i{n,i}); %i_max
        extreme_point{n,i}(2,1)=min(mult_ao_j{n,i}); %j_min
        extreme_point{n,i}(2,2)=max(mult_ao_j{n,i}); %j_max
    end
end

%%
for n=1:num2
    for i=1:abs(bubble{n,2})
        mult_flag{n,i}(1,1)=0; %是否优弧
        mult_flag{n,i}(1,2)=0; %类型
        mult_flag{n,i}(2,1)=0; %是否拟合
        mult_flag{n,i}(3,1)=0; %连线对应弧段
    end
    for i=1:abs(bubble{n,2})-1
        %判断是否为优弧,(1,1)=1表示该弧是优弧
        if(length(mult_ao_i{n,i})>major_arc_th)
            OA=[double(ao_point_de2{n,i}(1,1))-mult_mid_point{n,i}(2,1),double(ao_point_de2{n,i}(1,2))-mult_mid_point{n,i}(2,2)];
            OB=[double(ao_point_de2{n,i+1}(1,1))-mult_mid_point{n,i}(2,1),double(ao_point_de2{n,i+1}(1,2))-mult_mid_point{n,i}(2,2)];
            if(acos(dot(OA,OB)/(norm(OA)*norm(OB)))/pi*180<major_arc_angle)
                mult_flag{n,i}(1,1)=1;
            else
                mult_flag{n,i}(1,1)=0;
            end
        end
        %判断弧的类型(1、凸向左下；2、右下；3、左上；4、右上)
        if(length(mult_ao_i{n,i})>major_arc_th)
            i_mult_line_mid_point=(mult_ao_i{n,i}(1,1)+mult_ao_i{n,i}(1,end))/2;
            j_mult_line_mid_point=(mult_ao_j{n,i}(1,1)+mult_ao_j{n,i}(1,end))/2;
            if(mult_mid_point{n,i}(2,1)>=i_mult_line_mid_point && mult_mid_point{n,i}(2,2)<=j_mult_line_mid_point)
                mult_flag{n,i}(1,2)=1;
            elseif(mult_mid_point{n,i}(2,1)>=i_mult_line_mid_point && mult_mid_point{n,i}(2,2)>=j_mult_line_mid_point)
                mult_flag{n,i}(1,2)=2;
            elseif(mult_mid_point{n,i}(2,1)<=i_mult_line_mid_point && mult_mid_point{n,i}(2,2)<=j_mult_line_mid_point)
                mult_flag{n,i}(1,2)=3;
            elseif(mult_mid_point{n,i}(2,1)<=i_mult_line_mid_point && mult_mid_point{n,i}(2,2)>=j_mult_line_mid_point)
                mult_flag{n,i}(1,2)=4;
            end
        end
    end
    i=abs(bubble{n,2});
    if(length(mult_ao_i{n,i})>major_arc_th)
        dis_A_mid(n,i)=norm([double(ao_point_de2{n,i}(1,1)),double(ao_point_de2{n,i}(1,2))]-[mult_mid_point{n,i}(2,1),mult_mid_point{n,i}(2,2)]);
        dis_B_mid(n,i)=norm([double(ao_point_de2{n,1}(1,1)),double(ao_point_de2{n,1}(1,2))]-[mult_mid_point{n,i}(2,1),mult_mid_point{n,i}(2,2)]);
        dis_A_B(n,i)=norm([double(ao_point_de2{n,i}(1,1)),double(ao_point_de2{n,i}(1,2))]-[double(ao_point_de2{n,1}(1,1)),double(ao_point_de2{n,1}(1,2))]);
        if(dis_A_mid(n,i)>=dis_A_B(n,i)&&dis_B_mid(n,i)>=dis_A_B(n,i))
            mult_flag{n,i}(1,1)=1;
        else
            mult_flag{n,i}(1,1)=0;
        end
    end
    if(length(mult_ao_i{n,i})>major_arc_th)
        i_mult_line_mid_point=(mult_ao_i{n,i}(1,1)+mult_ao_i{n,i}(1,end))/2;
        j_mult_line_mid_point=(mult_ao_j{n,i}(1,1)+mult_ao_j{n,i}(1,end))/2;
        if(mult_mid_point{n,i}(2,1)>i_mult_line_mid_point && mult_mid_point{n,i}(2,2)<j_mult_line_mid_point)
            mult_flag{n,i}(1,2)=1;
        elseif(mult_mid_point{n,i}(2,1)>i_mult_line_mid_point && mult_mid_point{n,i}(2,2)>j_mult_line_mid_point)
            mult_flag{n,i}(1,2)=2;
        elseif(mult_mid_point{n,i}(2,1)<i_mult_line_mid_point && mult_mid_point{n,i}(2,2)<j_mult_line_mid_point)
            mult_flag{n,i}(1,2)=3;
        elseif(mult_mid_point{n,i}(2,1)<i_mult_line_mid_point && mult_mid_point{n,i}(2,2)>j_mult_line_mid_point)
            mult_flag{n,i}(1,2)=4;
        end
    end
end
%%
%记录所有弧段对应的弧线序号
partner=zeros(num2,mult);
for n=1:num2
    for i=1:abs(bubble{n,2})
        i_mult_line_mid_point=(mult_ao_i{n,i}(1,1)+mult_ao_i{n,i}(1,end))/2;
        j_mult_line_mid_point=(mult_ao_j{n,i}(1,1)+mult_ao_j{n,i}(1,end))/2;
        i1=mult_mid_point{n,i}(2,1);
        j1=mult_mid_point{n,i}(2,2);
        i2=i_mult_line_mid_point;
        j2=j_mult_line_mid_point;
        e=1:5:1000;
        phi=atan((i1-i2)/(j1-j2));
        if(mult_flag{n,i}(1,2)==1)
            iiiii=i1+e.*sin(phi);
            jjjjj=j1+e.*cos(phi);
        elseif(mult_flag{n,i}(1,2)==2)
            iiiii=i1-e.*sin(phi);
            jjjjj=j1-e.*cos(phi);
        elseif(mult_flag{n,i}(1,2)==3)
            iiiii=i1+e.*sin(phi);
            jjjjj=j1+e.*cos(phi);
        elseif(mult_flag{n,i}(1,2)==4)
            iiiii=i1-e.*sin(phi);
            jjjjj=j1-e.*cos(phi);
        else
            continue;
        end
        num3=bubble{n,2};
        [~,num4]=size(e);
        for q=1:num4
            for w=1:num3
                if(i==w)
                    continue;
                elseif(iiiii(1,q)>extreme_point{n,w}(1,1)&&iiiii(1,q)<extreme_point{n,w}(1,2)&&...
                        jjjjj(1,q)>extreme_point{n,w}(2,1)&&jjjjj(1,q)<extreme_point{n,w}(2,2))
                partner(n,i)=w;
                break;
                end
            end
            if(partner(n,i)~=0)
                break;
            end
        end
    end
end
%%
%计算每条弧的平均清晰度
Qabf_result = Qabf(image);
mult_ao_Qabf=zeros(num2,mult);
for n=1:num2
    if(abs(bubble{n,2})>2)
        for i=1:abs(bubble{n,2})
            sum_arc=0;
            [~,num_mult]=size(mult_ao_i{n,i});
            for q=1:num_mult
                sum_arc=sum_arc+Qabf_result(mult_ao_i{n,i}(1,q),mult_ao_j{n,i}(1,q));
            end
            mult_ao_Qabf(n,i)=sum_arc/num_mult;
        end
    end
end
%%
%先尝试与partner拟合
for n=1:num2
    for q=1:(abs(bubble{n,2})-1)
        if(mult_flag{n,q}(2,1)~=1&&length(mult_ao_i{n,q})>15)
            if(partner(n,q)~=0)
                temp_i=[mult_ao_i{n,q}(1,:),mult_ao_i{n,partner(n,q)}(1,:)];
                temp_j=[mult_ao_j{n,q}(1,:),mult_ao_j{n,partner(n,q)}(1,:)];
                [temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,5)]=my_ellipsefit(temp_j,temp_i);
                cluster_result=ifcluster(temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,5),temp_i,temp_j,partner_cluster_distance,partner_cluster_proportion);
                if(cluster_result && abs(mult_ao_Qabf(n,q)-mult_ao_Qabf(n,w))/mult_ao_Qabf(n,q)<0.2)
                    mult_flag{n,q}(2,1) = 1;
                    mult_flag{n,partner(n,q)}(2,1) = 1;
                    bubble{n,6}{mult_finish(n,1),1}(1,:) = temp_ellipse(1,:);
                    bubble{n,9}{mult_finish(n,1),1}(:,1) = temp_i;    % 2021/06/28 添加
                    bubble{n,9}{mult_finish(n,1),1}(:,2) = temp_j;    % 2021/06/28 添加
                    mult_finish(n,1) = mult_finish(n,1)+1;
                    continue;
                end
            end
        end
    end
end

%对优弧进行椭圆拟合
for n=1:num2
    if(abs(bubble{n,2})>2)
        for i=1:abs(bubble{n,2})
            if(isempty(mult_flag{n,i})~=1 && mult_flag{n,i}(1,1)==1 )
                [mult_flag{n,i}(4,1), mult_flag{n,i}(4,2), mult_flag{n,i}(4,3), ...
                    mult_flag{n,i}(4,4), mult_flag{n,i}(4,5)] =  ...
                    my_ellipsefit(mult_ao_j{n,i}(1,1:length(mult_ao_j{n,i})),mult_ao_i{n,i}(1,1:length(mult_ao_i{n,i})));
                mult_flag{n,i}(2,1)=1;
                bubble{n,6}{mult_finish(n,1),1}=[mult_flag{n,i}(4,1), mult_flag{n,i}(4,2), ...
                    mult_flag{n,i}(4,3),mult_flag{n,i}(4,4), mult_flag{n,i}(4,5)];
                bubble{n,9}{mult_finish(n,1),1}(:,1) = mult_ao_i{n,i}(1,1:length(mult_ao_i{n,i}));    % 2021/06/28 添加
                bubble{n,9}{mult_finish(n,1),1}(:,2) = mult_ao_j{n,i}(1,1:length(mult_ao_j{n,i}));    % 2021/06/28 添加
                mult_finish(n,1)=mult_finish(n,1)+1;
            end
        end
    end
end

%优弧拟合出的椭圆与其余圆弧聚类
for n=1:num2
    if (abs(bubble{n,2})>2 & isempty(bubble{n,6}))
        [b,~]=size(bubble{n,6});
        for q=1:b
            if(mult_flag{n,q}(2,1)==1)
                for w=1:abs(bubble{n,2})
                    if(abs(bubble{n,2})<4)
                        if(abs(mult_ao_Qabf(n,q)-mult_ao_Qabf(n,w))/mult_ao_Qabf(n,q)<0.2&&ifcluster(mult_flag{n,q}(4,3),mult_flag{n,q}(4,4),mult_flag{n,q}(4,1), ...
                                mult_flag{n,q}(4,2),mult_flag{n,q}(4,5),mult_ao_i{n,w},mult_ao_j{n,w},0.05,0.95))
                        mult_flag{n,w}(2,1)=1;
                        end
                    else
                        if(abs(mult_ao_Qabf(n,q)-mult_ao_Qabf(n,w))/mult_ao_Qabf(n,q)<0.2&&ifcluster(mult_flag{n,q}(4,3),mult_flag{n,q}(4,4),mult_flag{n,q}(4,1), ...
                                mult_flag{n,q}(4,2),mult_flag{n,q}(4,5),mult_ao_i{n,w},mult_ao_j{n,w},0.2,0.9))
                        mult_flag{n,w}(2,1)=1;
                        end
                    end
                end
            end
        end
    end
end
%%
%剩下的两两聚合
temp_ellipse=zeros(1,5);
temp_i=zeros();
temp_j=zeros();
for n=1:num2
    if(abs(bubble{n,2})>2)
        for q=1:(abs(bubble{n,2})-1)
            if(mult_flag{n,q}(2,1)~=1&&length(mult_ao_i{n,q})>10)
                %先尝试与partner拟合
                if(partner(n,q)~=0)
                    if(mult_flag{n,partner(n,q)}(2,1)~=1)
                        temp_i=[mult_ao_i{n,q}(1,:),mult_ao_i{n,partner(n,q)}(1,:)];
                        temp_j=[mult_ao_j{n,q}(1,:),mult_ao_j{n,partner(n,q)}(1,:)];
                        [temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,5)]=my_ellipsefit(temp_j,temp_i);
                        if(abs(bubble{n,2})<4)
                            cluster_result=ifcluster(temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,5),temp_i,temp_j,0.05,0.95);
                        else
                            cluster_result=ifcluster(temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,5),temp_i,temp_j,0.2,0.8);
                        end
                        if(cluster_result&&abs(mult_ao_Qabf(n,q)-mult_ao_Qabf(n,w))/mult_ao_Qabf(n,q)<0.2)
                            mult_flag{n,q}(2,1)=1;
                            mult_flag{n,partner(n,q)}(2,1)=1;
                            bubble{n,6}{mult_finish(n,1),1}(1,:)=temp_ellipse(1,:);
                            bubble{n,9}{mult_finish(n,1),1}(:,1) = temp_i;    % 2021/06/28 添加
                            bubble{n,9}{mult_finish(n,1),1}(:,2) = temp_j;    % 2021/06/28 添加
                            mult_finish(n,1)=mult_finish(n,1)+1;
                            continue;
                        end
                    end
                end
                %再两两聚合
                for w=(q+1):abs(bubble{n,2})
                    if(length(mult_ao_i{n,w})>10 &&mult_flag{n,w}(2,1)~=1&&mult_flag{n,q}(1,2)~=mult_flag{n,w}(1,2))
                        temp_i=[mult_ao_i{n,q}(1,:),mult_ao_i{n,w}(1,:)];
                        temp_j=[mult_ao_j{n,q}(1,:),mult_ao_j{n,w}(1,:)];
                        [temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,5)]=my_ellipsefit(temp_j,temp_i);
                        if(abs(bubble{n,2})<4)  
                            cluster_result=ifcluster(temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,5),temp_i,temp_j,0.05,0.95);
                        else
                            cluster_result=ifcluster(temp_ellipse(1,3),temp_ellipse(1,4),temp_ellipse(1,1),temp_ellipse(1,2),temp_ellipse(1,5),temp_i,temp_j,0.2,0.8);
                        end
                        if(cluster_result&&abs(mult_ao_Qabf(n,q)-mult_ao_Qabf(n,w))/mult_ao_Qabf(n,q)<0.2)
                            mult_flag{n,q}(2,1)=1;
                            mult_flag{n,w}(2,1)=1;
                            bubble{n,6}{mult_finish(n,1),1}(1,:)=temp_ellipse(1,:);
                            bubble{n,9}{mult_finish(n,1),1}(:,1) = temp_i;    % 2021/06/28 添加
                            bubble{n,9}{mult_finish(n,1),1}(:,2) = temp_j;    % 2021/06/28 添加
                            mult_finish(n,1)=mult_finish(n,1)+1;
                        end
                    end
                    if(mult_flag{n,q}(2,1)==1)
                        continue;
                    end
                end
            end
        end
    end
end

%剩余拟合
for n=1:num2
    for q=1:abs(bubble{n,2})
        if(mult_flag{n,q}(2,1)~=1 && length(mult_ao_i{n,q})>=50)
            [bubble{n,6}{mult_finish(n,1),1}(1,1),bubble{n,6}{mult_finish(n,1),1}(1,2),bubble{n,6}{mult_finish(n,1),1}(1,3), ...
                bubble{n,6}{mult_finish(n,1),1}(1,4),bubble{n,6}{mult_finish(n,1),1}(1,5)]=my_ellipsefit(mult_ao_j{n,q}(1,:),mult_ao_i{n,q}(1,:));
            bubble{n,9}{mult_finish(n,1),1}(:,1) = mult_ao_i{n,q}(1,:);    % 2021/06/28 添加
            bubble{n,9}{mult_finish(n,1),1}(:,2) = mult_ao_j{n,q}(1,:);    % 2021/06/28 添加
            mult_finish(n,1)=mult_finish(n,1)+1;
        end
    end
end

%剩余拟合  不附加条件，所有弧段一起拟合
for n=1:num2
    re_i = zeros();
    re_j = zeros();
    for q=1:abs(bubble{n,2})
        if(mult_flag{n,q}(2,1)~=1)
            re_i = [re_i,mult_ao_i{n,q}(1,:)];
            re_j = [re_j,mult_ao_j{n,q}(1,:)];
        end
    end
    re_i(:,1)=[];
    re_j(:,1)=[];
    [bubble{n,6}{mult_finish(n,1),1}(1,1),bubble{n,6}{mult_finish(n,1),1}(1,2),bubble{n,6}{mult_finish(n,1),1}(1,3), ...
    bubble{n,6}{mult_finish(n,1),1}(1,4),bubble{n,6}{mult_finish(n,1),1}(1,5)]=my_ellipsefit(re_j,re_i);
    bubble{n,9}{mult_finish(n,1),1}(:,1) = re_i;    % 2021/06/28 添加
    bubble{n,9}{mult_finish(n,1),1}(:,2) = re_j;    % 2021/06/28 添加
    mult_finish(n,1)=mult_finish(n,1)+1;
end

% output:
bub_overlap = bubble;

end



function [result] = Qabf(image)
% 输出图片清晰度
% model parameters
% L=1; Tg=0.9994;kg=-15;Dg=0.5;Ta=0.9879;ka=-22;Da=0.8;       
% Sobel Operator
h1=[1 2 1;0 0 0;-1 -2 -1]; h2=[-1 0 1;-2 0 2;-1 0 1];
% if y is the response to h1 and x is the response to h3;
% then the intensity is sqrt(x^2+y^2) and orientation is arctan(y/x);
if size(image,3)==3
    image=rgb2gray(image);
end
image = im2double(image);
SAx = conv2(image,h2,'same'); SAy = conv2(image,h1,'same');
gA = sqrt(SAx.^2 + SAy.^2); 
% [M,N] = size(SAx); aA = zeros(M,N);
result=gA;
end

function [ifclusterresult] = ifcluster(x0,y0,a,b,phi,i,j,paramater1,paramater2)

if(nargin<8)
    paramater1=0.2;
    paramater2=0.8;
end
c=sqrt(a^2-b^2);
x1=x0+c*cos(phi);
y1=y0-c*sin(phi);
x2=x0-c*cos(phi);
y2=y0+c*sin(phi);
num_cluster_point=0;
for n=1:length(i)
    dist=sqrt((j(1,n)-x1)^2+(i(1,n)-y1)^2)+sqrt((j(1,n)-x2)^2+(i(1,n)-y2)^2);
    if(abs(dist-2*a)/(2*a)<paramater1)
        num_cluster_point=num_cluster_point+1;
    end
end
if(real(a)<10||real(b)<10)
    ifclusterresult=0;
elseif(num_cluster_point/length(i)>paramater2)
    ifclusterresult=1;
else
    ifclusterresult=0;
end
end