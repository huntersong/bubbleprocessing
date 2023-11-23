%参数定义
%图像坐标(x,y)和像素坐标(u,v)，世界坐标(U,V,Z)，_world坐标系，_camera坐标系  _info全信息
image_camera1 = imread("figure1.png"); %假定每个图像中气泡占用的像素小于128px；
image_camera1 = rgb2gray(image_camera1);
image_camera1 = im2double(image_camera1);
[y,x] = size(image_camera1);
image_camera2 = zeros(y,x);
U = x; V = y; W = max(x,y); 
bubble_face_position = zeros(U*V*W,3);  %气泡空间坐标（x，y，z）
bubble_center_world = [0,0,0];  %气泡空间中心坐标（x，y，z）
bubble_voxel_store = zeros(U*V*W,3);
miu1 = 1;  miu2 = 1.3;     %折射率介质1  折射率介质2
rho_1_2 = miu1/miu2;   
para_world = [0,0,0];  %世界坐标系
para_camera = [0,0,0];  %相机参数
u_camera1_px = 128;  %相机图像坐标(x,y)和像素坐标(u,v)
v_camera1_px = 128;
U_camera1_px = 128;  %相机世界坐标(x,y)和像素坐标(U,V)
V_camera1_px = 128;
Z_camera1_px = 128;
U_LED1_face_world = 1;  %LED的世界坐标(x,y)和像素坐标(U,V)
V_LED1_face_world = 2;
Z_LED1_face_world = 3;
theta_face1_world = 1;
phi_face1_world =2;
u_camera2_px_camera = 128;
v_camera2_px_camera = 128;
x_LED2_face_world = 1;
y_LED2_face_world = 2;
z_LED1_face_world = 3;
x_world_bubble = 0;
y_world_bubble = 0;
z_world_bubble = 0;
%图例信息
ray_position = [0,0,0];
ray_direction = [0,0,0];
ray_info_world = [1,ray_position,ray_direction];  % 
ray_position_image = [0,0,0];
ray_direction_image = [0,0,0];
d_imagetobubble = 0;
ray_image_world = [1,ray_position_image,ray_direction_image,d_imagetobubble];%   c点  i1矢量  d1：image到气泡前面的距离
bubble_face_position = [0,0,0]; 
ray_face_direction = [0,0,0]; bubble_ray_info_world = [1,ray_face_direction];    %   v1点   o1  
theta_face_world = 0; phi_face_world = 0; 
bubble_face_direction = [theta_face_world,phi_face_world];   bubble_faceinfo_world = [0,0,0];     %   n1
bubble_position_face_world = [bubble_face_position,bubble_faceinfo_world];
ray_position_LED = [0,0,0];
ray_direction_LED = [0,0,0];
d_bubbletoLED = 0;
ray_LED_world = [1,ray_position_LED,ray_direction_LED,d_bubbletoLED];   %   q点   o2矢量  d2：LED到气泡后面的距离
omega_camera = 10/pi; %相机弥散元内视角锥
omega_LED = 10/pi; %漫射光源光锥

%光源参数：光源和一般任意光源不同，

%程序

%（1）{标定过程
n = 10;  m = 10;    %标定板中特征点的数目
position_calibration = zeros(n,m);    %标定板中特征点的坐标

for n=1
    step=1;
    %标定：1、找到标定板中每个点的图像的像素坐标（u_c,v_c)
    % 读入左右两张图像
    left_img = imread('left.png');
    right_img = imread('right.png');

    % 棋盘格参数
    num_corners = [8, 6];
    square_size = 20; % 棋盘格边长

    % 棋盘格角点检测
    [imagePoints, boardSize] = detectCheckerboardPoints({left_img, right_img}, 'NumRows', num_corners(1), 'NumColumns', num_corners(2));

    % 生成世界坐标系下的棋盘格点云
    worldPoints = generateCheckerboardPoints(boardSize, square_size);

    % 单目标定
    [params_left, ~, ~] = estimateCameraParameters(imagePoints{1}, worldPoints);
    [params_right, ~, ~] = estimateCameraParameters(imagePoints{2}, worldPoints);

    % 双目标定
    [rotationMatrix, translationVector, ~, ~, ~, reprojectionErrors] = estimateWorldCameraPose(imagePoints{1}, worldPoints, params_left);
    stereoParams = stereoParameters(params_left, params_right, rotationMatrix, translationVector);

    % 保存结果
    save('stereoParams.mat', 'stereoParams');
    % 标定的结果存储为：
    % 建立的坐标系->以相机1为原点
    camera1_world = [0,0,0];%原点   
    camera2_world = [100,0,-500];    camera3_world = [0,0,0];
    LED1_world = [0,0,-600];    LED2_world = [-80,0,-200];
    P_a_camera1and2 = [25,147,-800]; P_b_camera1and2 = [238,147,-1000] ;...
        P_c_camera1and2 = [238,147,-1400]; P_d_camera1and2 = [25,147,-1050];
end    

%}
%（2）{

    %视体算法,双相机钝角视体算法，silhouette method，钝角可以拍摄部分气泡背侧的信息，

    gray_img_camera1 = image_camera1;     gray_img_camera2 = image_camera2;   %灰度图
    canny_img_camera1 = edge(gray_img_camera1, 'canny');    canny_img_camera2 = edge(gray_img_camera2, 'canny');    %二值化和边缘空洞填充
    imfill_img_camera1 = imfill(canny_img_camera1, 'holes'); imfill_img_camera2 = imfill(canny_img_camera2, 'holes');
    %point1_img_camera1 = transform2dto3d(1,1);  point2_img_camera1 = transform2dto3d(1,1);% 根据标定结果将（x，y）转换到（x，y，z）
    point1_img_camera1 = [25,147,-1000];  point2_img_camera1 = [238,147,-1100];  ...
        point1_img_camera2 = [100,147,-900];  point2_img_camera2 = [80,147,-1200];      %每一层气泡切片的四个顶点 三个相机是六个顶点
    
    bubble_center_world = 1/4 * (point1_img_camera1 + point2_img_camera1 + point1_img_camera2 + point2_img_camera2);    %计算每一层气泡的中心点
    
    % silhouette方法 把顶点用弧线连接起来
    bubble_up = 16;     bubble_down = 278;    
    for n=1:1    % 这是第147层的气泡切片
        [xx_selected_1,yy_selected_1] = bezier_points(point1_img_camera1,P_a_camera1and2,point1_img_camera2);
        [xx_selected_2,yy_selected_2] = bezier_points(point1_img_camera2,P_b_camera1and2,point2_img_camera1);
        [xx_selected_3,yy_selected_3] = bezier_points(point2_img_camera2,P_c_camera1and2,point2_img_camera1);
        [xx_selected_4,yy_selected_4] = bezier_points(point1_img_camera1,P_d_camera1and2,point2_img_camera2);
        zz_selected_1 = n*ones(1,length(xx_selected_1));  zz_selected_2 = n*ones(1,length(xx_selected_2));...
            zz_selected_3 = n*ones(1,length(xx_selected_3));  zz_selected_4 = n*ones(1,length(xx_selected_4));
        bubble_silhouette = [xx_selected_1',yy_selected_1',zz_selected_1';xx_selected_2',yy_selected_2',zz_selected_2';...
            xx_selected_3',yy_selected_3',zz_selected_3';xx_selected_4',yy_selected_4',zz_selected_4'];   % 这是以相机1为原点，建立的坐标系
        plot(bubble_silhouette(:,1), bubble_silhouette(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);axis equal
    end

    %bubble_silhouette = [count,U,V,Z,1]; 
    ptCloud = pointCloud(bubble_silhouette);   ptCloud;     %从输入点坐标创建一个点云对象。检查点云对象的属性
    pcshow(ptCloud);
    normals = pcnormals(ptCloud);   ptCloud = pointCloud(bubble_silhouette,'Normal',normals); %在每个点上增加面法向：pcnormals评估这个点云上每个点的法向量
    pcshow(ptCloud);
    U_bubble_face_world = ptCloud.Location(:,1);    V_bubble_face_world = ptCloud.Location(:,2);...
        W_bubble_face_world = ptCloud.Location(:,3);  %点坐标
    U_normals_bubble_face_world = normals(:,1);    V_normals_bubble_face_world = normals(:,2);...
        W_normals_bubble_face_world = normals(:,3); %构建气泡表面的法向量
    hold on;    
    quiver3(U_bubble_face_world,V_bubble_face_world,W_bubble_face_world,U_normals_bubble_face_world,V_normals_bubble_face_world,W_normals_bubble_face_world);    
    hold off

    %（约束气泡，空间点（ xi_w，yi_w，zi_w ）， ）
%}
%（3）
%基于光追的细观结构{
%光路初级判断，
image_camera1 = imread("figure1.png"); image_01_camera1 = rgb2gray(double(image_camera1));
slice_left = 25; slice_right = 238; slice_up = 1;  slice_down = 1;  %获得气泡二值化图像的上下左右每一层的边缘，以减小光追的计算量
for v_camera_px = 1:128     % 相机 1 图像上 气泡所占水平方向像素的空间若为100，则增加到128来覆盖这个重建区域
    for u_camera_px = 1:128     %索寻光路贯通（灰度值预判）
%find Imax的u，v坐标 （u，v为最左上角点，分四个象限）
        [M,I] = max(gray_img_camera1(147,:)); % 找到矩阵 A 气泡的每一层中的最大值和位置
        [I_row, I_col] = ind2sub(size(gray_img_camera1), I); 
        gray_max_image1 = M;
    end
end
maxValue_x_cam1 = I_row;   maxValue_y_cam1 = I_col;
bubble_equivalent_diameter = 200/2;  %先暂定一个等价直径（需要再考虑）
d_maxValue_imagetobubble = [bubble_center_world(1,1), bubble_center_world(1,2),...
    bubble_center_world(1,3) + bubble_equivalent_diameter];  %bubble_equivalent_diameter这里是给定的假定值，后面将调整

ray_count_max = 100;%中心高光点信息：射线数n（=100；I=100k）
weight = gray_max_image1/255;%中心高光点信息：
intensity_ray_max = ray_count_max*weight;%中心高光点信息：射线数n（=100；I=100k）
theta_face_world = 0; phi_face_world = 0;   %中心高光点面法向信息：气泡面偏转角度
x_bubble_info_max_theta = 0; y_bubble_info_max_theta = 0 ; z_bubble_info_max_theta = 1;%中心高光点面法向信息：气泡面偏转角度
bubble_faceinfo_world = zeros(length(bubble_silhouette),3);  %构建中心高光点信息：气泡前面偏转角度 矩阵
ray_face_direction = zeros(length(bubble_silhouette),3);    %气泡面光线折射信息o1
bubble_position_face_world = [bubble_silhouette,bubble_faceinfo_world]; 

bubble_position_face_world(maxValue_x_cam1-slice_left+1,4:6) = [x_bubble_info_max_theta, ...
    y_bubble_info_max_theta, z_bubble_info_max_theta]; %中心高光点信息：最亮点的光路前面信息
U_maxValue_world = bubble_position_face_world(maxValue_x_cam1-slice_left+1,1); V_maxValue_world = bubble_position_face_world(maxValue_x_cam1-slice_left+1,1);
%在当前切片上的二维角度通过 z_bubble_info_max_theta来确定 ；x_bubble_info_max_theta = 0, y_bubble_info_max_theta = 0
%初始光坐标 i1 这里光强定了100*weight，
%分四个象限进行索引，首先是第1象限,按maxValue_x_cam1 = u_camera_px;   maxValue_y_cam2 = v_camera_px
%区域Ⅰ开始索引，随后Ⅱ，Ⅲ，Ⅳ区域  结合silhouette结果将灰度图像首次关联到三维气泡的表面
for v_camera_px = 1:(slice_up-maxValue_y_cam1+1)  %区域Ⅰ开始索引
    for u_camera_px = 1:(slice_right-maxValue_x_cam1)             
        U_bubble_info_world = (u_camera_px - 1) + U_maxValue_world; 
        V_bubble_info_world = (v_camera_px - 1) + V_maxValue_world; 
        W_bubble_info_world = bubble_center_world(1,3); 
        z_theta = pi/2*(1 - gray_img_camera1(u_camera_px,v_camera_px)/gray_max_image1);     % 二维情况下的表面倾角
        temp_bubble_info_world = [U_bubble_info_world, V_bubble_info_world, W_bubble_info_world, x_theta = 0, y_theta = 0 , z_theta]  %结合silhouette结果将灰度图像首次关联到三维气泡的表面
        %射线数/灰度值：nI=100，nI=60（基于几何光学的光线追踪分析） 
    %气液界面法向量与射线夹角：θ=0（ nI=100 ），则表示界面与相机面平行，
    %例如：判断？当θ=30°还是-30°（ nI=60和θ ）。   
    % 构建三维情况下（θ，φ） 
        theta_face_world = z_theta;    %和z轴夹角
        phi_face_world = arcsin(x_theta/sqrt(x_theta^2+y_theta^2));     %x-y平面和y夹角 在（4）中进一步修正
        bubble_info_world(m,:) = temp_bubble_info_world;    %气液界面的世界坐标信息
        m = m + 1;  
    end
end
%}

%（4）
% {
%分四个象限进行索引，首先是第1象限,按maxValue_x_cam1 = u_camera_px;   maxValue_y_cam2 = v_camera_px
%区域Ⅰ开始索引，随后Ⅱ，Ⅲ，Ⅳ区域             %按照ppt19页进行推进空间索引。
for region = 1:4
    region = 1; % 从区域Ⅰ开始
    bubble_faceinfo_world = bubble_info_world;
    thetatemp = theta_face_world;
    m = 1;
    for u_camera_px = 1:(128-ray_ori_px(1,1))
        for v_camera_px = 1:(128-ray_ori_px(1,2))
            if thetatemp(m)>26.57
                %下个空间_world索引区域为a
                bubble_faceinfo_world(m,2) = bubble_faceinfo_world(m,2)+1;  % 修正 U_bubble_info_world, V_bubble_info_world, W_bubble_info_world； 修正phi_face_world
                bubble_faceinfo_world(m,3) = bubble_faceinfo_world(m,3)-1;
            elseif -26.57<thetatemp<26.57
                    %下个索引区域为b
                bubble_faceinfo_world(m,2) = bubble_faceinfo_world(m,2)+1;  % 修正 U_bubble_info_world, V_bubble_info_world, W_bubble_info_world； 修正phi_face_world
                bubble_faceinfo_world(m,3) = bubble_faceinfo_world(m,3);
            else %下个索引区域为c
                %A网格 α=0.5   B网格 θ=30°和交点，α=1-tanθ
                %已知气泡（ xi_w，yi_w，zi_w ），界面起点A（ xi，yi，zi ）
                %		则B（ xi+1，yi-1，zi ）
                bubble_faceinfo_world(m,2) = bubble_faceinfo_world(m,2)+1;  % 修正 U_bubble_info_world, V_bubble_info_world, W_bubble_info_world； 修正phi_face_world
                bubble_faceinfo_world(m,3) = bubble_faceinfo_world(m,3)+1;               
            end
            %存储为（ xi+1，yi-1，zi ），α，（θ，φ）
                    %->
            % {  %光追程序：寻找真正的前面法向n1矢量：n1矢量的优化过程 ->  o1矢量 v1
            %c点光线的数目=射线数n；
            [bubble_info_world, theta_face_world, phi_face_world] = raytrace_adjust(ray_image_world , bubble_faceinfo_world, ray_LED_world);
            %对该点的信息进行存储`
            theta_face_world = ;    %气液界面的法向角度
            phi_face_world = ;      %
            bubble_info_world = [bubble_face_position,ray_face_direction];   %气液界面的世界坐标
            ray_info = [ray_position,ray_direction];    %气泡前界面处的光路信息
        end

    end    
    % }
end

% }综上为一个气泡正面光追，对于背侧的气泡重建，将正面与背面对换重新迭代。








