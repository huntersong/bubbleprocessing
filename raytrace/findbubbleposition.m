function findbubbleposition
%把外面的白色变成黑色
image_camera1 = imread('D:\办公电脑文件\reynoldsstress\Motion Detection\raytrace\raytracebuble\figure1.png');
if numel(size(image_camera1))>2
    disp('  彩色度图像\n');
    image_camera1 = rgb2gray(image_camera1);
end
image_camera1 = im2double(image_camera1);
[y,x] = size(image_camera1);
image_camera1 = imcrop(image_camera1,[2 2 (x-2*2) (y-2*2)]);

[m,n]=size(imcrop_bubble);
for x=1:m
    for y=1:n
        if(imcrop_bubble(x,y)==255)
            imcrop_bubble(x,y)=0;
        end
    end
end
filename = [sprintf('%02d',n),'data.mat'];
save(filename,bubble_all_position_face_world);



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
bub_file_path = my_read('D:\办公电脑文件\reynoldsstress\Motion Detection\raytrace\raytrace','png');
image_back = my_background('D:\办公电脑文件\reynoldsstress\Motion Detection\raytrace\raytrace\background_00001.png'); 
image_back = uint8(255/65535*double(image_back));
i = 2;
image_bub = imread(bub_file_path{i});
image_bub = uint8(255/65535*double(image_bub));
[image,bub_bw,gray_image,lab_image] = my_pre_processing(image_bub,image_back,'imagesize',[1,1500,1,1000]); 
location = my_sobel(bub_bw,im2double(gray_image),lab_image);

end

