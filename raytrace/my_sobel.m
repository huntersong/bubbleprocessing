function out = my_sobel(int_bw,int_im,int_L,varargin)
% 本函数的作用是为气泡图像处理算法提供失焦气泡去除功能，使用方式如下：
% 
% 输入：
%     int_bw：原始气泡边界信息
%     int_im：原始气泡图像
%     int_L：带标记的气泡位置图
%     varargin：可变长度输入列表，用于指定阈值阈值
% 
% 输出：
%     out：输出为失焦气泡的位置信息location
% 
% 调用说明：
%     out = my_sobel(int_bw,int_im,int_L); 
%                  默认预处理方式，输出为失焦气泡位置信息
%     out = my_sobel(int_bw,int_im,int_L,'bubble-gray-th',0.35,'bubble-size-th',50);
%                  指定气泡灰度阈值为0.35，阈值范围为[0,1];指定气泡边界梯度判定尺寸阈值为50，尺寸越大，可判定气泡越多，可能会造成错误；
% 
% 版本号VOL1.0，编写于2021年6月3日，作者：WG-Chen
% 模糊处理、识别

% default value
bub_gray_th = 0.25;
bub_grad_size_th = 50;
grad_th = 0.25;

if numel(varargin) > 0
    for i = 1 : numel(varargin)
        tf1 = strcmp(varargin{i},'bubble-grad-th');
        tf2 = strcmp(varargin{i},'bubble-size-th');
        if tf1 == 1
            grad_th = str2double(varargin{i+1});
        end
        if tf2 == 1
            bub_grad_size_th = str2double(varargin{i+1});
        end
    end
end

% origrayimage = im2double(int_im);
origrayimage = int_im;
[im_S,im_L] = size(int_im);

[objectnum,~] = size(int_bw);
bou_nei_gray = cell(objectnum,2);

for i = 1:objectnum
    bounlength = size(int_bw{i,1});
    for j = 1:bounlength
         for ii = 1:3
            bou_nei_gray{i,1}(3*j-2,ii) = origrayimage(int_bw{i,1}(j,1)-1,int_bw{i,1}(j,2)+ii-2);
            bou_nei_gray{i,1}(3*j-1,ii) = origrayimage(int_bw{i,1}(j,1),int_bw{i,1}(j,2)+ii-2);
            bou_nei_gray{i,1}(3*j,ii) = origrayimage(int_bw{i,1}(j,1)+1,int_bw{i,1}(j,2)+ii-2);
         end
    end
end

%定义sobel算子
sobelx = [-1,-2,-1;0,0,0;1,2,1];
sobely = sobelx';

%开始计算边界的梯度
for i = 1 :objectnum
    bounlength = size(int_bw{i,1});
    for j = 1:bounlength
     %定义移动窗口 slidwin
     slidwin = [bou_nei_gray{i,1}(3*j-2,1),bou_nei_gray{i,1}(3*j-2,2),bou_nei_gray{i,1}(3*j-2,3);...
         bou_nei_gray{i,1}(3*j-1,1),bou_nei_gray{i,1}(3*j-1,2),bou_nei_gray{i,1}(3*j-1,3);bou_nei_gray{i,1}(3*j,1),...
         bou_nei_gray{i,1}(3*j,2),bou_nei_gray{i,1}(3*j,3)];
     %转换为double型
     slidwin = double(slidwin);
     Sx = sobelx .* slidwin;
     Sy = sobely .* slidwin;
     %注意，这里为了提高效率，可以使用绝对值算数平均数来代替几何平均数
     %      gradx = abs(sum(sum(Sx)));
     %      grady = abs(sum(sum(Sy)));
     %      bou_nei_gray{i,2}(j,1) = gradx + grady;
     %准确的方法应当算x、y方向的几何平均数
     gradx = sum(sum(Sx));
     grady = sum(sum(Sy));
     bou_nei_gray{i,2}(j,1) = sqrt(gradx^2+grady^2);
    end
end
%计算曲线灰度梯度最大值、最小值，并记录在数组中
maxgrad = zeros(objectnum,1);
mingrad = zeros(objectnum,1);
meangrad = zeros(objectnum,1);
bub_gray = zeros(objectnum,3); % 1 mean 2 area_num 3 max

for i = 1 : im_S
    for j = 1 : im_L
        if int_L(i,j) ~= 0
            bub_gray(int_L(i,j),1) = bub_gray(int_L(i,j),1) + origrayimage(i,j);
            bub_gray(int_L(i,j),2) = bub_gray(int_L(i,j),2) + 1 ;
            if bub_gray(int_L(i,j),3) == 0
                bub_gray(int_L(i,j),3) = origrayimage(i,j);
            elseif origrayimage(i,j) < bub_gray(int_L(i,j),3)
                bub_gray(int_L(i,j),3) = origrayimage(i,j);
            end
        end
    end
end

for i = 1:objectnum
    maxgrad(i,1) = max(bou_nei_gray{i,2});
    mingrad(i,1) = min(bou_nei_gray{i,2});
    meangrad(i,1) = mean(bou_nei_gray{i,2});
    bub_gray(i,1) = bub_gray(i,1)/bub_gray(i,2);
end

%计算最大梯度
gromaxgrad = max(meangrad);
bub_grad_th = gromaxgrad*grad_th;
graymax = max(bub_gray(:,1));
gray_th = graymax*0.25;

location = zeros();
grad_num = 0;
[A,~] = size(meangrad);

for g_i =  1 : A
    if meangrad(g_i,1) < bub_grad_th & length(int_bw{g_i,1}) < bub_grad_size_th
        grad_num = grad_num + 1 ;
        location(grad_num,1) = g_i;
    elseif bub_gray(g_i,3) > bub_gray_th
        grad_num = grad_num + 1 ;
        location(grad_num,1) = g_i;
%     elseif bub_gray(g_i,1) < gray_th
%         grad_num = grad_num + 1 ;
%         location(grad_num,1) = g_i;
    end
end

out = location;
end
% 
% image_1 = int_im;
% for i = 1 : length(bubble_boundaries)
%     for j = 1 : length(bubble_boundaries{i,1})
%         image_1(bubble_boundaries{i,1}(j,1),bubble_boundaries{i,1}(j,2)) = 255;
%     end
% end
% image_1 = int_im;
% for i = 1 : length(int_bw)
%     for j = 1 : length(int_bw{i,1})
%         image_1(int_bw{i,1}(j,1),int_bw{i,1}(j,2)) = 255;
%     end
% end
% figure
% imshow(image_1);
% hold on
% for i = 1 : length(location)
%     text(int_bw{location(i,1),1}(1,2),int_bw{location(i,1),1}(1,1),[num2str(location(i,1)),'(',num2str(meangrad(i,1)),')'],'color','b');
% end
% for i = 1 : length(int_bw)
%     if meangrad(i,1) < grad_th & length(int_bw{i,1}) < 40
%         text(int_bw{i,1}(1,2),int_bw{i,1}(1,1),[num2str(i),'(',num2str(meangrad(i,1)),')'],'color','b');
%     else
%         text(int_bw{i,1}(1,2),int_bw{i,1}(1,1),[num2str(i),'(',num2str(meangrad(i,1)),')'],'color','r');
%     end
% end
% hold off

