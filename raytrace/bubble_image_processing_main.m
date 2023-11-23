% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 本软件为气泡图像提供数字处理功能：
% 平均计算速度：238秒/百张，约2.4秒/张
% 重要说明：
% bubble内存储了识别到的气泡信息：1 是否斑点；2 凹点数量；3.圆度； 4.面积； 5.长宽比； 6.气泡拟合信息；7.拟合气泡所用的边界；
% bubble_overlap内存储了重叠气泡的信息： 1 是否斑点；2 凹点数量；3.圆度； 4.面积； 5.长宽比； 6.气泡拟合信息；
%                                                7.该重叠边界在bubble上的位置 8.该重叠边界所有的分段弧 9.拟合气泡所用的边界；
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 版本号VOL1.0，编写于2021年6月5日，作者：WG-Chen
% 修改日期 2021/06/28 增加内容：优化重复气泡的删除判定；增加识别凹点后分段弧；增加重叠气泡所用的识别弧段信息；
% ！！重要提示：结果显示，算法中取凹点后分段弧存在误差，主要集中在起点和终点的选择
% ！！重要提示：结果显示，可能是拟合的方式还没有得到优化，待优化
% 修改日期 2021/07/08 增加内容：my_bubble_tracking.m 用于做气泡的二维、三维追踪；
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close;clear;clc;

% % 指定存储路径
% path = uigetdir(cd,'请选择气泡识别结果存储路径');
% 获取气泡图像存储路径   
bub_file_path = my_read('D:\办公电脑文件\reynoldsstress\Motion Detection\raytrace\raytrace','png');  
% 'F:\算法\bubble-image-processing\基于机器学习&分段弧聚类的重叠气泡处理方案','bmp'
%'E:\连续谱气泡发生器实验\0601CWG\1200sccm,300sccm,8.2278\2000~high','bmp'
% 'F:\课题组PPT\0608图像处理算法PPT\基于高速相机的一站式多相流测试方案\微流控\bmp','bmp'
% 'G:\bmp\New folder','bmp'
% 识别背景图    
image_back = my_background('D:\办公电脑文件\reynoldsstress\Motion Detection\raytrace\raytrace\background_00001.png');   
%'F:\算法\bubble-image-processing\基于机器学习&分段弧聚类的重叠气泡处理方案\figure\test_zhang_back.bmp'
%'E:\连续谱气泡发生器实验\0601CWG\background\back_00375.bmp'
% 'C:\Users\Admin\Desktop\yxy20210709\back1.bmp'

% h = waitbar(0,'正在计算中，请稍等！');

for i = 1 : numel(bub_file_path)
    
    clearvars -except i bub_file_path image_back cutting_point h bubble_tracking_data
    
    image_bub = imread(bub_file_path{i});
%     if i == 1
%         [image,bub_bw,gray_image,lab_image,cutting_point] = my_pre_processing(image_bub,image_back);   %,'imagesize',[1,1018,1,2029]
%     else
%         [image,bub_bw,gray_image,lab_image,~] = my_pre_processing(image_bub,image_back,'cutting',cutting_point);
%     end

    [image,bub_bw,gray_image,lab_image] = my_pre_processing(image_bub,image_back,'imagesize',[1,1500,1,1000]);   %[500,1530,1,2000][70,860,350,750]
    
    location = my_sobel(bub_bw,im2double(gray_image),lab_image);
    if location~=0
        bub_bw([location(:,1)],:)=[];
    end
%     % 进度条
%     waitbar(i/numel(bub_file_path),h,['正在计算中，请稍等！(',num2str(i),'/',num2str(numel(bub_file_path)),')']);
    
    [bubble,bub_overlap,overlap_ao_data] = my_bub_processing(bub_bw,lab_image);
    
    if size(bub_overlap,1) ~= 0
        bubble_overlap = my_overlapbubbles(bub_overlap,overlap_ao_data,image);

        [num,~] = size(bubble_overlap);
        for j =  1 : num
            for n = 1 : numel(bubble_overlap{j,6})
                bubble{bubble_overlap{j,7},6}{n,1} =  bubble_overlap{j,6}{n,1};      % 重叠气泡识别结果
                bubble{bubble_overlap{j,7},7}{n,1} =  bubble_overlap{j,9}{n,1};      % 重叠气泡识别所用的弧
            end
        end
    end
    
    bubble = my_bub_deleting(bubble,'maxsize_th',1000);   %
    
    if i == 1
        bubble_tracking_data = my_bubble_tracking(bubble,'first');
    else
        bubble_tracking_data = my_bubble_tracking(bubble,'others',bubble_tracking_data,'maxlength_th',20);
    end
    
% 用于实时作图，每张耗时约1s，可注释
    my_plot(image,bubble); 
    path_name_bubble = [bub_file_path{i,1}(1:length(bub_file_path{i,1})-4),'result.bmp'];
%     saveas(gcf,path_name_bubble);
%     close Figure 1
    
    bubble_data = my_postprocessing(bubble);
    
    bubble_size = bubble_data(:,2);
    bubble_out = my_statistics(bubble_size);
    close Figure 1
    
%     path_name_data = [bub_file_path{i,1}(1:length(bub_file_path{i,1})-4),'.xlsx'];
%     xlswrite(path_name_data,bubble_data);

end



figure;   % 气泡轨迹
imshow(image);
hold on
n = 5;
for i = 1 : bubble_tracking_data{n,1}
    plot_ellipse(bubble_tracking_data{n,3}{1,i}{1,3}(1,1),bubble_tracking_data{n,3}{1,i}{1,3}(1,2),...
        bubble_tracking_data{n,3}{1,i}{1,4}(1,1),bubble_tracking_data{n,3}{1,i}{1,4}(1,2),bubble_tracking_data{n,3}{1,i}{1,4}(1,3));
    text(bubble_tracking_data{n,3}{1,i}{1,3}(1,1),bubble_tracking_data{n,3}{1,i}{1,3}(1,2),'*','fontsize',20,'color','r');
    hold on
    if i  == 1
        pause(10)
    else
        pause(1)
    end
end


% figure;   %  画带标记的气泡
% imshow(image);
% hold on;
% [num2,~] = size(bubble);
% for n=1:num2
%     if(bubble{n,2}>1)
%         [b,~]=size(bubble{n,6});
%         for i=1:b
%             if(bubble{n,6}{i,1}~=0)
%                 plot_ellipse(real(bubble{n,6}{i,1}(1,3)),real(bubble{n,6}{i,1}(1,4)),real(bubble{n,6}{i,1}(1,1)),real(bubble{n,6}{i,1}(1,2)),real(bubble{n,6}{i,1}(1,5)));
%                 text(real(bubble{n,6}{i,1}(1,3)),real(bubble{n,6}{i,1}(1,4)),[num2str(n),'(',num2str(i),')']);
%                 hold on;
%             end
%         end
%     elseif(bubble{n,2}==1 & bubble{n,6}~=0)
%         plot_ellipse(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,2)),real(bubble{n,6}(1,5)));
%         text(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),num2str(n));
%         hold on;
%     elseif(bubble{n,2}==0.5)
%         plot_ellipse(real(bubble{n,6}(2,2)),real(bubble{n,6}(2,1)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,1)),0);
%         text(real(bubble{n,6}(2,2)),real(bubble{n,6}(2,1)),num2str(n));
%         hold on;
%     elseif(bubble{n,2}==0 & bubble{n,6}~=0)
%         plot_ellipse(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),real(bubble{n,6}(1,1)),real(bubble{n,6}(1,2)),real(bubble{n,6}(1,5)));
%         text(real(bubble{n,6}(1,3)),real(bubble{n,6}(1,4)),num2str(n));
%         hold on;
%     end
% end

% for n = 1 : numel(bubble{12,7})   % 某气泡拟合所用的弧
% plot(bubble{12,7}{n,1}(:,2),bubble{12,7}{n,1}(:,1),'linewidth',3);
% pause(3)
% end

% for n = 1 : numel(bubble_overlap{4,8}{1,1})    % 重叠边界凹点分割后的分段弧
% 
%     plot(bubble_overlap{4,8}{1,1}{n,1}(:,2),bubble_overlap{4,8}{1,1}{n,1}(:,1),'linewidth',2);
%     pause(3)
% end

% aa = imread('F:\算法\bubble-image-processing\基于机器学习&分段弧聚类的重叠气泡处理方案\test_zhang.bmp');
% aa = rgb2gray(aa);
% aa = im2double(aa);
% imwrite(aa,'test_zhang.bmp');
% aaa = ones(701,903);
% imwrite(aaa,'test_zhang_back.bmp');

% for i = 1 : 17   % rgb -> gray
%     aa = imread(['F:\课题组PPT\0608图像处理算法PPT\基于高速相机的一站式多相流测试方案\微流控\西交大 微流控 5KF10M 1280×860@1000fps_Moment',num2str(i),'.jpg']);
%     aa = rgb2gray(aa);
%     aa = im2double(aa);
%     imwrite(aa,['F:\课题组PPT\0608图像处理算法PPT\基于高速相机的一站式多相流测试方案\微流控\bmp\西交大 微流控 5KF10M 1280×860@1000fps_Moment',num2str(i),'.bmp']);
%     clear aa
% end
