function bubble_data = my_postprocessing_GAN(bubble,varargin)
% 本函数的作用是为气泡图像处理算法提供气泡后处理，使用方式如下：
% 
% 输入：
%     bubble：重叠气泡边界信息
%     N×6的元胞数组，{1}存储是否为斑点气泡；{2}存储凹点个数；{3}为圆度；{4}为面积；{5}为长宽比；{6}储存为拟合信息
%     varargin：可变长度输入列表，用于指定阈值
% 
% 输出：
%     bubble_data：后处理输出的气泡识别结果
% 
% 调用说明：
%     bubble_data = my_postprocessing(bubble)
%                  默认后处理方式，输出为识别的气泡结果
%     bubble_data = my_postprocessing(bubble,varargin)
%                  设定圆度为0.94
% 
% 版本号VOL1.0，编写于2021年6月4日，作者：WG-Chen

%% default value



%% methods


%% 计算部
[num2,~] = size(bubble);
d_bubble = struct;
k = 0;
for n=1:num2
    if(abs(bubble{n,2})==0.5 & bubble{n,6}~=0)
        k=k+1;
        d_bubble(k).aa = bubble{n,6}(1,1);
        d_bubble(k).bb = bubble{n,6}(1,1);
        d_bubble(k).phi = 0;
        d_bubble(k).xx = bubble{n,6}(2,1);
        d_bubble(k).yy = bubble{n,6}(2,2);
        d_bubble(k).BubDiameter = 2*bubble{n,6}(1,1);
        d_bubble(k).BubBoundary = bubble{n,7};

    elseif(abs(bubble{n,2})==0 & bubble{n,6}~=0)
        k=k+1;
        d_bubble(k).aa = bubble{n,6}(1,1);
        d_bubble(k).bb = bubble{n,6}(1,2);
        d_bubble(k).phi = bubble{n,6}(1,5);
        d_bubble(k).xx = bubble{n,6}(1,4);
        d_bubble(k).yy = bubble{n,6}(1,3);
        d_bubble(k).BubDiameter = ((2*abs(bubble{n,6}(1,1)))^2*2*abs(bubble{n,6}(1,2)))^(1/3);
        d_bubble(k).BubBoundary = bubble{n,7};

    elseif(abs(bubble{n,2})==1 & bubble{n,6}~=0)
        k=k+1;
        d_bubble(k).aa = bubble{n,6}(1,1);
        d_bubble(k).bb = bubble{n,6}(1,2);
        d_bubble(k).phi = bubble{n,6}(1,5);
        d_bubble(k).xx = bubble{n,6}(1,4);
        d_bubble(k).yy = bubble{n,6}(1,3);
        d_bubble(k).BubDiameter = ((2*abs(bubble{n,6}(1,1)))^2*2*abs(bubble{n,6}(1,2)))^(1/3);
        d_bubble(k).BubBoundary = bubble{n,7};

    elseif(abs(bubble{n,2})==2)
        if(real(bubble{n,6}{1,1}(1,1))~=0 && real(bubble{n,6}{1,1}(1,2))~=0)
            k=k+1;
            d_bubble(k).aa = bubble{n,6}{1,1}(1,1);
            d_bubble(k).bb = bubble{n,6}{1,1}(1,2);
            d_bubble(k).phi = bubble{n,6}{1,1}(1,5);
            d_bubble(k).xx = bubble{n,6}{1,1}(1,4);
            d_bubble(k).yy = bubble{n,6}{1,1}(1,3);
            d_bubble(k).BubDiameter = ((2*abs(bubble{n,6}{1,1}(1,1)))^2*2*abs(bubble{n,6}{1,1}(1,2)))^(1/3);
            d_bubble(k).BubBoundary = bubble{n,7}{1,1};
        end
        if(real(bubble{n,6}{2,1}(1,1))~=0 && real(bubble{n,6}{2,1}(1,2))~=0)
            k=k+1;
            d_bubble(k).aa = bubble{n,6}{2,1}(1,1);
            d_bubble(k).bb = bubble{n,6}{2,1}(1,2);
            d_bubble(k).phi = bubble{n,6}{2,1}(1,5);
            d_bubble(k).xx = bubble{n,6}{2,1}(1,4);
            d_bubble(k).yy = bubble{n,6}{2,1}(1,3);
            d_bubble(k).BubDiameter = ((2*abs(bubble{n,6}{2,1}(1,1)))^2*2*abs(bubble{n,6}{2,1}(1,2)))^(1/3);
            d_bubble(k).BubBoundary = bubble{n,7}{2,1};
        end
    elseif(abs(bubble{n,2})>=2)
        [b,~]=size(bubble{n,6});
        for i=1:b
            if(isempty(bubble{n,6}{i,1})~=1 & bubble{n,6}{i,1}~=0)
                k=k+1;
                d_bubble(k).aa = bubble{n,6}{i,1}(1,1);
                d_bubble(k).bb = bubble{n,6}{i,1}(1,2);
                d_bubble(k).phi = bubble{n,6}{i,1}(1,5);
                d_bubble(k).xx = bubble{n,6}{i,1}(1,4);
                d_bubble(k).yy = bubble{n,6}{i,1}(1,3);
                d_bubble(k).BubDiameter = ((2*abs(bubble{n,6}{i,1}(1,1)))^2*2*abs(bubble{n,6}{i,1}(1,2)))^(1/3);
                d_bubble(k).BubBoundary = bubble{n,7}{i,1};
            end
        end
    end
end

% outputs:
bubble_data = d_bubble;

end