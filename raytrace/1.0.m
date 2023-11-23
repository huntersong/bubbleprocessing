

for n = 147
    region = 1; % 从区域Ⅰ开始
    nn = n-1+point1_img_camera1(1,2);
    stats = regionprops('table', image_01_camera1(nn,:), 'BoundingBox','Extrema');
    boundbox = stats.BoundingBox;
    slice_left = round(boundbox(1,1));  slice_right = slice_left+round(boundbox(1,3));
    slice_up = 1;  slice_down = 1;  %获得气泡这一层的二值化图像的上下左右每一层的边缘，以减小光追的计算量
    [M,I] = max(gray_img_camera1(nn,:)); % 找到矩阵 A 气泡的每一层中的最大值和位置
    [I_row, I_col] = ind2sub(size(gray_img_camera1(nn,:)), I);
    gray_max_image1 = max(max(gray_img_camera1(nn,:)));
    maxValue_x_cam1 = I_col;
    maxValue_y_cam1 = I_row;

    ray_count_max = 100;    %中心高光点信息：射线数n（=100；I=100k）
    weight = gray_max_image1/255;   %中心高光点信息：
    intensity_ray_max = ray_count_max*weight;   %中心高光点信息：射线数n（=100；I=100k）
    theta_face_world = 0; phi_face_world = 0;   %中心高光点面法向信息：气泡面偏转角度
    x_bubble_info_max_theta = 0; y_bubble_info_max_theta = 0 ; z_bubble_info_max_theta = 1;%中心高光点面法向信息：气泡面偏转角度

    %   进行bubble_faceinfo_world 面参数修正;
    bubbletemp_silhouette = bubble_all_position_face_world(tempab(n,1):tempab(n,2),1:3);%气泡正面初始
    bubbletemp_back_silhouette = bubble_all_position_face_world((tempab(n,2)+1):(tempab(n,2)-tempab(n,1)+tempab(n,2)-1),1:3);%气泡背面初始
    bubbletemp_face = bubble_all_face_direction(tempab(n,1):tempab(n,2),1:3);%气泡正面法向

    bubble_faceinfo_world = zeros(length(bubbletemp_silhouette),3);  %构建中心高光点信息：气泡前面偏转角度 矩阵
    bubble_face_direction = zeros(length(bubbletemp_silhouette),3);  %气泡面法向
    ray_face_direction = zeros(length(bubbletemp_silhouette),3);    %气泡面光线折射信息o1
    bubble_position_face_world = [bubbletemp_silhouette,bubble_faceinfo_world];

    bubble_position_face_world(maxValue_x_cam1-slice_left+1,4:6) = [x_bubble_info_max_theta, ...
        y_bubble_info_max_theta, z_bubble_info_max_theta]; %中心高光点信息：最亮点的光路前面信息
    U_maxValue_world = bubble_position_face_world(maxValue_x_cam1-slice_left+1,1); ...
        V_maxValue_world = bubble_position_face_world(maxValue_x_cam1-slice_left+1,2);...
        W_maxValue_world = bubble_position_face_world(maxValue_x_cam1-slice_left+1,3);

    thetatemp = bubble_all_position_face_world(tempab(n,1):tempab(n,2),6);  %temp_z_theta
    m = 1;
    for v_camera_px = 1:(slice_up-maxValue_y_cam1+1)
        for u_camera_px = 1:(slice_right-maxValue_x_cam1)
            U_bubble_info_world = U_maxValue_world + (u_camera_px - 1); 
            V_bubble_info_world = V_maxValue_world + (v_camera_px - 1); 
            threshold = 10;  %value 应根据灰度梯度差来不停的修正
            temp_z_theta = pi/2*(1 - gray_img_camera1(V_bubble_info_world,U_bubble_info_world)/gray_max_image1);
            value = threshold*temp_z_theta;
            m = U_bubble_info_world-slice_left+1;
            if thetatemp(m)>(26.57*pi/180)
                %下个空间_world索引区域为a    
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3)-value;% 修正W_bubble_info_world坐标位置； 修正phi_face_world
            elseif thetatemp(m)<(26.57*pi/180) && thetatemp(m)>-(26.57*pi/180)
                %下个索引区域为b 
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3); % 修正W_bubble_info_world坐标位置； 修正phi_face_world   
            else %下个索引区域为c
                %A网格 α=0.5   B网格 θ=30°和交点，α=1-tanθ
                %已知气泡（ xi_w，yi_w，zi_w ），界面起点A（ xi，yi，zi ）
                %		则B（ xi+1，yi-1，zi ）
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3)+value;% 修正W_bubble_info_world坐标位置； 修正phi_face_world               
            end
            %存储为（ xi+1，yi-1，zi ），α，（θ，φ）
                    %->
            % {  %光追程序：寻找真正的前面法向n1矢量：n1矢量的优化过程 ->  o1矢量 v1
            %c点光线的数目=射线数n；
            input_bubbletemp_silhouette = bubbletemp_silhouette(m,1:3);
            input_bubbletemp_face = bubbletemp_face(m,1:3);

            [out_bubbletemp_silhouette,output_bubbletemp_face] = ...
                raytrace_adjust(input_bubbletemp_silhouette,input_bubbletemp_face);
            value1 = 0;  %后续考虑一个修正值的确定方法
            bubbletemp_silhouette(m,3) = out_bubbletemp_silhouette(m,3)+value1;
            %对该点的信息进行存储
            bubble_position_face_world(m,1:3) = bubbletemp_silhouette(m,1:3);            
        end
        for u_camera_px = 1:(maxValue_x_cam1-slice_left)
            U_bubble_info_world = U_maxValue_world - (u_camera_px-1);
            V_bubble_info_world = V_maxValue_world + (v_camera_px - 1);
            threshold = 10;  %value 应根据灰度梯度差来不停的修正
            temp_z_theta = pi/2*(1 - gray_img_camera1(V_bubble_info_world,U_bubble_info_world)/gray_max_image1);
            value = threshold*temp_z_theta;
            m = U_bubble_info_world-slice_left+1;
            if thetatemp(m)>(26.57*pi/180)
                %下个空间_world索引区域为a    
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3)-value;% 修正W_bubble_info_world坐标位置； 修正phi_face_world
            elseif thetatemp(m)<(26.57*pi/180) && thetatemp(m)>-(26.57*pi/180)
                %下个索引区域为b 
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3); % 修正W_bubble_info_world坐标位置； 修正phi_face_world   
            else %下个索引区域为c
                %A网格 α=0.5   B网格 θ=30°和交点，α=1-tanθ
                %已知气泡（ xi_w，yi_w，zi_w ），界面起点A（ xi，yi，zi ）
                %		则B（ xi+1，yi-1，zi ）
                bubbletemp_silhouette(m,3) = bubbletemp_silhouette(m,3)+value;% 修正W_bubble_info_world坐标位置； 修正phi_face_world               
            end
            %存储为（ xi+1，yi-1，zi ），α，（θ，φ）
                    %->
            % {  %光追程序：寻找真正的前面法向n1矢量：n1矢量的优化过程 ->  o1矢量 v1
            %c点光线的数目=射线数n；
            input_bubbletemp_silhouette = bubbletemp_silhouette(m,1:3);
            input_bubbletemp_face = bubbletemp_face(m,1:3);

            [out_bubbletemp_silhouette,output_bubbletemp_face] = ...
                raytrace_adjust(input_bubbletemp_silhouette,input_bubbletemp_face);
            value1 = 0;  %后续考虑一个修正值的确定方法
            bubbletemp_silhouette(m,3) = out_bubbletemp_silhouette(m,3)+value1;
            %对该点的信息进行存储
            bubble_position_face_world(m,1:3) = bubbletemp_silhouette(m,1:3); 
        end
    end
    bubble_all_position_face_world(tempab(n,1):tempab(n,2),1:3) = bubble_position_face_world(:,1:3);
end

plot3(bubble_all_position_face_world(44929:45145,1),bubble_all_position_face_world(44929:45145,2),...
    bubble_all_position_face_world(44929:45145,3),'k-o')
hold on
plot3(bubble_all_position_face_world(45146:45362,1),bubble_all_position_face_world(45146:45362,2),...
    bubble_all_position_face_world(45146:45362,3),'k-o')
axis equal
view(0,180)



t=3;            %LED的面积    
sparsity = 2;   %可以理解为LED内灯泡密度，或者是散射点密度，sparsity数值越小越密集
planeNormal_LED = [0 0 1];
planePoint_LED = [0 0 -1500];
ray_position_LED = planePoint_LED(1,3)*ones(x*y*t^2,3);   %
mm=1;
for yy = (0-y*(t-1)/2):sparsity:(0+y*(t+1)/2)
    for xx = (0-y*(t-1)/2):sparsity:(0+y*(t+1)/2)
        ray_position_LED(mm,1) = xx;
        ray_position_LED(mm,2) = yy;
        mm = mm+1;
    end
end

ray_direction_LED = [0 0 1];
    d_bubbletoLED = ;
    ray_LED_world = [ray_position_LED, ray_direction_LED, d_bubbletoLED]; %q点  o2向量







