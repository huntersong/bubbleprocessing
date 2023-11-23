function image_back = my_background(path)
% 本函数的作用是为气泡图像处理算法提供读取气泡背景图，使用方式如下：
% 
% 输入：
%     path：背景图路径，可不输入
% 
% 输出：
%     image_back：输出背景图
% 
% 调用说明：
%     image_back = my_background(); 
%                  无路径输入，跳出文件选择对话框
%     image_back = my_background(path)
%                  根据输入路径读取背景图像文件
% 
%版本号VOL1.0，编写于2021年6月2日，作者：WG-Chen

if nargin == 0
    [background_file,background_path] = uigetfile('*.*','请选择背景：');
    background_path = fullfile(background_path,background_file);
    image_back = imread(background_path);
end
    
if nargin == 1
    image_back = imread(path);
end    

end