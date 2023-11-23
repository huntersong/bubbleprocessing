function [imfill_img_camera1] = binary_tempfunc(image)

image_camera1 = image; %假定每个图像中气泡占用的像素小于128px；
if numel(size(image_camera1))>2
    fprintf('%s','彩色度图像');
    image_camera1 = rgb2gray(image_camera1);
end
image_camera1 = im2double(image_camera1);
[y,x] = size(image_camera1);
%image_camera1 = imcrop(image_camera1,[2 2 (x-2*2) (y-2*2)]);

I_gray = image_camera1;%原图像变为灰度图像
level=graythresh(I_gray);%计算图像I_gray的全局阈值，level为标准化灰度值，其范围为[0 1]
[height,width]=size(I_gray);%计算灰度图像的长宽
I_bw=im2bw(I_gray,0.7);

for i=1:height %%循环中进行反色
    for j=1:width   
        if I_bw(i,j)==1      
            I_bw(i,j)=0;  
        else I_bw(i,j)=1; 
        end
    end
end
%figure(1);imshow(I_bw);%显示取反后的二值图像（背景为黑色）
I_bw = imfill(I_bw,'holes');
I_double=double(I_bw);
I_double=edge(I_double,'sobel');
imfill_img_camera1 = I_bw;
imfill_img_camera2 = I_bw;

end



