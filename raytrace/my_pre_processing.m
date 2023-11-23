function [out_image,bub_bw,ori_image,lab_image,cutting_point] = my_pre_processing(image_bub,image_back,varargin)
% 本函数的作用是为气泡图像处理算法提供图像预处理功能，使用方式如下：
% 
% 输入：
%     image_bub：气泡图像
%     image_back：背景图像
%     varargin：可变长度输入列表，用于选择预处理方法、阈值
% 
% 输出：
%     image：去背景后的气泡图像
%     bub_bw：识别到的气泡边界信息
%     ori_image：切割后的气泡原图，用于删除失焦气泡，详见my_sobel
%     lab_image：带标签的气泡位置图，用于删除失焦气泡，详见my_sobel
%     cutting_point：切割点信息
% 
% 调用说明：
%     bub_bw = my_pre_processing(image_bub,image_back); 
%                  默认预处理方式，输出为识别的边界信息
%     [image,bub_bw,ori_image,lab_image] = my_pre_processing(image_bub,image_back,...
%                                         'smallsize',30,'Gaussian_Filtering','imagesize',[1,1000,1,1000],'cutting'); 
%                  过滤小泡尺寸设定为30pixels，选择高斯双核滤波，设定输出增加原图和带标签位置图
% 
% 版本号VOL1.0，编写于2021年6月2日，作者：WG-Chen

%% default value
small_bub_size = 20;
Gaussian_Filtering = 0;
flag = 0;
flag_back = 1;

%% methods
if image_back == 0
    flag_back = 0;
end


if numel(varargin) == 0
    small_bub_size = 20;
    Gaussian_Filtering = 0;
else
    for i = 1 : numel(varargin)
        tf1 = strcmp(varargin{i},'smallsize');
        tf2 = strcmp(varargin{i},'Gaussian_Filtering');
        tf3 = strcmp(varargin{i},'imagesize');
        tf4 = strcmp(varargin{i},'cutting');
        if tf1 == 1
            small_bub_size = varargin{i+1};
        end
        if tf2 == 1
            Gaussian_Filtering = 1;
        end
        if tf3 == 1 
            cutting_point = varargin{i+1};
            flag = 1;
        end
        if tf4 == 1
            get_point = varargin{i+1};
            flag = 2;
        end
    end
end
    
if flag == 0
    % 取点
    get_point = my_getpoint(image_bub);
    get_point = round(get_point);
    image_bub_cut = my_cutting(get_point,image_bub);
    if flag_back
        image_back_cut = my_cutting(get_point,image_back);
    end
elseif flag == 1
    image_bub_cut = image_bub(cutting_point(1):cutting_point(2),cutting_point(3):cutting_point(4));
    if flag_back
        image_back_cut = image_back(cutting_point(1):cutting_point(2),cutting_point(3):cutting_point(4));
    end
else 
    image_bub_cut = my_cutting(get_point,image_bub);
    if flag_back
        image_back_cut = my_cutting(get_point,image_back);
    end
end

%% 计算部
ori_image = image_bub_cut;

bwbubblemax = max(max(image_bub_cut));
image_bub_fix = image_bub_cut.*((2^8)/bwbubblemax); 
if flag_back
    image_back_fix = image_back_cut.*((2^8)/bwbubblemax); 
    image = image_back_fix - image_bub_fix;
    out_image = image;
else
    out_image = image_bub_fix;
    image = image_bub_fix;
end


image = imadjust(image);                               % 灰度调整
if Gaussian_Filtering == 1                             % 高斯双核滤波，去噪保边
    image = imbilatfilt(image);
end
image = medfilt2(image,[3,3]);                         % 中值滤波
gy = graythresh(image);                                % 大津算法确定二值化阈值
image = imbinarize(image,gy);                          % 二值化处理
% image = 1 - image;
image_fill = imfill(image,'holes');                    % 填洞处理
image_fill = bwareaopen(image_fill,small_bub_size);    % 去除小泡
[image_label,~] = bwlabel(image_fill,8);               % 对气泡区域图标号
image_label = my_del_boundbub(image_label);            % 删除与边界粘连的气泡
bubble_boundaries = bwboundaries(image_label);         % 取气泡边界信息

for i = 1 : numel(bubble_boundaries)                   % 将原始气泡边界编号记录在bubble_boundaries{i,2}中
    bubble_boundaries{i,2} = i;
end

bub_bw = bubble_boundaries;
lab_image = image_label;

if flag ~= 1 
    cutting_point = get_point;
end

end

% plot画边界 plot(bubble_boundaries{1,1}(:,2),bubble_boundaries{1,1}(:,1))


%% 以下为本函数中使用到的子函数

function point = my_getpoint(image)
% 本函数用于图像取点，取一个长方形四个角点，取点顺序为：左上、左下、右上、右下点

h = msgbox({'Please select the segmentation point in the following figure, in the order of: upper left, lower left, upper right and lower right point;',...
    'The method of picking points is: click a point in the image with the mouse and press enter to record the picking points.'},'Select a point in the diagram');
set(h,'Position',[500 500 500 150]);
% 修改字体
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 15); 

pause(2)

point = zeros(4,2);
image = adapthisteq(image);
figure
imshow(image)
set(gcf,'outerposition',get(0,'screensize'));
for j = 1 : 4
    [point(j,1),point(j,2)] = ginput;
end
close all

end


function image_label_out = my_del_boundbub(image_label_in)
% 本函数用于去除带标号连通域内与图像边界联通的区域
[row,column] = size(image_label_in);
L = image_label_in;
for i=1:row                                                               %去除与边界粘贴的气泡
    if(L(i,1)~=0)
        c=L(i,1);
        cc=find(L==c);
        L(cc)=0;
    end
end
for k=1:6
    for i=1:row
        if(L(i,column-k+1)~=0)
            c=L(i,column-k+1);
            cc=find(L==c);
            L(cc)=0;
        end
    end
end

for j=1:column
    if(L(1,j)~=0)
        c=L(1,j);
        cc=find(L==c);
        L(cc)=0;
    end
end
for j=1:column
    if(L(row,j)~=0)
        c=L(row,j);
        cc=find(L==c);
        L(cc)=0;
    end
end

[L2,~] = bwlabel(L,8);                                                 %去除与边界重合的气泡后重新对联通区域标号

image_label_out = L2;

end