function [position_x,position_y,position_z] = findbubbleposition(image_camera1,image_background)
%image_background = imread('background_00001.png');image_camera1 = imread('case_27406.png');
if numel(size(image_camera1))>2
    disp('  彩色度图像\n');
    image_camera1 = rgb2gray(image_camera1);
end
image_camera1 = double(image_camera1);
[y,x] = size(image_camera1);
image_camera1 = imcrop(image_camera1,[2 2 (x-2*2) (y-2*2)]);
if max(max(image_camera1))>255
    disp("image is > 8 bit")
    image_camera1 = uint8(255/65535*double(image_camera1));
else
    disp("image is 8 bit")
end 

if numel(size(image_background))>2
    disp('  彩色度图像\n');
    image_background = rgb2gray(image_background);
end
image_background = double(image_background);
[y,x] = size(image_background);
image_background = imcrop(image_background,[2 2 (x-2*2) (y-2*2)]);
if max(max(image_background))>255
    disp("image is > 8 bit")
    image_background = uint8(255/65535*double(image_background));
else
    disp("image is 8 bit")
end 

C=double(image_camera1)-double(image_background);   %image_background=imadjust(image_background,[0,1],[1,0]);
C=mat2gray(C);
Z2=imbinarize(C,0.1);
Z2=medfilt2(Z2,[4,4]);      %img2 = ~Z2;
I_gray = image_camera1;
level=graythresh(I_gray);%计算图像I_gray的全局阈值，level为标准化灰度值，其范围为[0 1]
[height,width]=size(I_gray);%计算灰度图像的长宽
I_bw=imbinarize(I_gray,level);%im2bw使用阈值level将灰度图像转换为二值图像．
for i=1:height %%循环中进行反色
    for j=1:width   
        if I_bw(i,j)==1      
            I_bw(i,j)=0;  
        else 
            I_bw(i,j)=1; 
        end
    end
end
% figure(3);imshow(I_bw);%显示取反后的二值图像（背景为黑色）

[L,num]=bwlabel(I_bw,8);%bwlabel标注二值图像I_bw中的目标物体，返回标识矩阵Ｌ和I_bw中目标物体的数量num，８表示连通数．
plot_x=zeros(1,num);%%zeros(m,n)产生m×n的全0矩阵．用于记录质心位置的横坐标
plot_y=zeros(1,num);%zeros(m,n)产生m×n的全0矩阵．用于记录质心位置的纵坐标

for k=1:num  %%num个区域依次统计质心位置 
    sum_x=0;    sum_y=0;    area=0; %初始化
    for i=1:height   
        for j=1:width 
            if L(i,j)==k     
                sum_x=sum_x+i;  %计算第Ｋ区域的横坐标总和
                sum_y=sum_y+j;  %计算第Ｋ区域的纵坐标总和 
                area=area+1;    %计算第Ｋ区域的由多少个坐标点表示
            end
        end
    end
    plot_x(k)=fix(sum_x/area);  %计算第Ｋ区域的质心横坐标
    plot_y(k)=fix(sum_y/area);%计算第Ｋ区域的质心纵坐标
end
position_x = plot_x(k);
position_y = plot_y(k);
position_z = 1;

[image,bub_bw,gray_image,lab_image] = my_pre_processing(image_camera1,image_background,'imagesize',[1,1497,1,997]); 

%{
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
bub_file_path = my_read('D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case5','bmp');
image_back = uint8(ones(301,351)*255);
i = 2;
image_bub = imread(bub_file_path{i});
[x y] = size(image_bub);
if max(max(image_bub))>255
    disp("image is > 8 bit")
    image_bub = uint8(255/65535*double(image_bub));
else
    disp("image is 8 bit")
end 

[image,bub_bw,gray_image,lab_image] = my_pre_processing(image_bub,image_back,'imagesize',[1,x,1,y]); 
location = my_sobel(bub_bw,im2double(gray_image),lab_image);
%}
end



