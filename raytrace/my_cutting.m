function image_out = my_cutting(input_point,input_image)
% 该函数用于将输出图像input_image根据输入的角点矩阵input_point切割，并输出相同大小的图像image_out
% 
% 输入：
%     input_point：切割图像的四个角点位置信息
%     input_image：带切割图像
% 
% 输出：
%     image_out：切割后的图像
% 
% 调用说明：
%     image_out = my_cutting(input_point,input_image);
%                  默认预处理方式，输出切割后的图像
% 
% 版本号VOL1.0，编写于2021年6月2日，作者：WG-Chen

image = input_image;
[L,S] = size(image);
ret_line = cell(4,1);
% im_line 函数用于获取两点之间的直线上的点集
ret_line{1,1} = my_line(input_point(1,1),input_point(1,2),input_point(2,1),input_point(2,2));
ret_line{2,1} = my_line(input_point(1,1),input_point(1,2),input_point(3,1),input_point(3,2));
ret_line{3,1} = my_line(input_point(3,1),input_point(3,2),input_point(4,1),input_point(4,2));
ret_line{4,1} = my_line(input_point(2,1),input_point(2,2),input_point(4,1),input_point(4,2));

% 下面将四个框框起来，抽取各标定板所占位置
% 首先将四条线连成一个完整的图像边界

% 边界连接顺序：左上->左下->右下->右上->左上

tran_bw = ret_line{1,1};
for j = 1 : length(ret_line{4,1})-1
    tran_bw(length(ret_line{1,1})+j,:) = ret_line{4,1}(j+1,:);
end
A = length(tran_bw);
B = 0;
for j = length(ret_line{3,1})-1 : -1 : 1
    B = B + 1;
    tran_bw(A+B,:) = ret_line{3,1}(j,:);
end
A = length(tran_bw);
B = 0;
for j = length(ret_line{2,1})-1 : -1 : 1
    B = B + 1;
    tran_bw(A+B,:) = ret_line{2,1}(j,:);
end
ret_bw = tran_bw;
clear A B tran_bw


% 抠图
tran_image = zeros(L,S);
for j = 1 : length(ret_bw)
    tran_image(ret_bw(j,1),ret_bw(j,2)) = 1;
end
tran_image = imfill(tran_image);
[A,B] = size(tran_image);
for n = 1 : A
    for m = 1 : B
        if tran_image(n,m) == 1
            image(n,m) = image(n,m);
        end
    end
end
clear tran_image A B


% output
image_out = image(min(input_point(:,2)):max(input_point(:,2)),min(input_point(:,1)):max(input_point(:,1)));

end



%% 以下为my_cutting函数的子函数，获取两个点连线中的所有点集
function Coor = my_line(x1,y1,x2,y2)
% 输入为x1,y1为点1，x2，y2为点2；输出Coor为识别到的两点之间直线上的所有点集
% x2-x1 不能等于零

if x2 == x1
    x2 = x2 + 1;
end

deltaH = abs(y1 - y2);
deltaW = abs(x1 - x2);
k = (y2-y1)/(x2-x1);
b = y1 - k*x1;
Num = 0;

if deltaH <= deltaW
    for j = x1 : x2
            Hb = y1 - 1;                            %以A点为起点
            He = y1 + 1;
            H = 0;
            W = 0;
            Min = 1000;        
            for i = Hb : He
                Tmpb = i - k*j;
                delta = abs(b - Tmpb);
                if delta < Min
                    Min = delta;
                    H = i;
                    W = j;
                end
            end
            if H ~= 0 && W ~= 0
                Num = Num + 1;                  %直线上像素点个数
                Coor(Num,1) = H;                %直线上像素点位置坐标
                Coor(Num,2) = W;
                y1 = H;
            end
    end
end

if deltaH > deltaW
    for i = y1 : y2
            Min = 1000;
            H = 0;
            W = 0;
            Wb = x1 - 4;                    %自变量变化范围可根据直线特征设定
            We = x1 + 4;
            for j = Wb : We
                Tmpb = i - k*j;
                delta = abs(b - Tmpb);
                if delta < Min
                    Min = delta;
                    H = i;
                    W = j;
                end
            end
            if H ~= 0 && W ~= 0
                Num = Num + 1;
                Coor(Num,1) = H;
                Coor(Num,2) = W;
                x1 = W;
            end
    end
end

end