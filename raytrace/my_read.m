function bub_file_path = my_read(path,image_type)
% 本函数的作用是为气泡图像处理算法提供读图功能、读取并记录图片路径，使用方式如下：
% 
% 输入：
%     path：气泡图片存储路径
%     image_type：气泡图像类型
% 
% 输出：
%     bub_file_path：气泡图像的完整路径名，1×N的元胞数组，N为检测到的气泡图像数
% 
% 调用说明：
%     bub_file_path = my_read(); 无路径输入，打开对话框选定气泡存储文件夹并输出检测到气泡路径
%     bub_file_path = my_read(path,image_type); 指定路径和气泡类型，输出气泡存储路径
% 
% 版本号VOL1.0，编写于2021年6月2日，作者：WG-Chen

if nargin == 2
    getfilename = ls([path,'\*.',image_type]');
    [num,~] = size(getfilename);
    if num == 0
        error('文件夹内未找到指定格式的气泡图像！');
    end
    filename = cellstr (getfilename);
    for i = 1 : num
        filename{i,1} = [path,'\',filename{i,1}];
    end
end

if nargin == 0
    path = uigetdir(cd,'请选择气泡图像路径：');
    liststring = {'bmp','png','jpg','tif','raw'};
    list = listdlg('PromptString','请选择气泡图像格式：','liststring',liststring,'ListSize',[150,100]);
    image_type = liststring{1,list};
    getfilename = ls([path,'\*.',image_type]');
    [num,~] = size(getfilename);
    if num == 0
        error('文件夹内未找到指定格式的气泡图像！');
    end
    filename = cellstr (getfilename);
    for i = 1 : num
        filename{i,1} = [path,'\',filename{i,1}];
    end
end


if nargout == 1
    bub_file_path = filename;
else
    error('输出变量数量只能为 1 ！'); 
end

end
