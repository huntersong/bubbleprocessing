function [output_bubbletemp_silhouette,output_bubbletemp_face] = raytrace_adjust(input_bubbletemp_silhouette,input_bubbletemp_face)    
%根据气泡图像和LED的信息来更新光线追踪的结果

planeNormal_LED = [0 0 1];
planePoint_LED = [0 0 -1500];

for ray_count = 1:100   %ray_number = 100;  
    miu1 = 1;  miu2 = 1.3;     %折射率介质1  折射率介质2
    rho_1_2 = miu1/miu2;
    %theta_ray_direction_image = (-5 + 10.*rand(1,1))*pi/180; %相机1的光锥；若在二维的角度，并设置光锥二维角度为10°，
    %phi_ray_direction_image = (-5 + 10.*rand(1,1))*pi/180;      %
        %ray_direction_image = [x_theta_image,y_theta_image,z_theta_image];
    z_theta_image = 1;
    x_theta_image = tan(pi/36)*(-1 + 2.*rand(1,1));
    y_theta_image = tan(pi/36)*(-1 + 2.*rand(1,1));
    ray_direction_LED = [x_theta_image,y_theta_image,z_theta_image];
    
    %立体角的国际制单位是球面度（steradian，sr）。立体角有一个非国际制单位平方度，1 sr = (180/π)2 
    
    %已知i1向量，c点， o2向量，q点（球坐标随机采样，追踪到LED？ ，与前文不同的是，
    %q点是不唯一的，导致o1也是不唯一的，这就决定了最终求得的n1法向量的不唯一性，
    %且n1的数目与光线相关，这也表示了CMOS的入射光线实际上也属于一个立体锥，
    %和LED的出射光一致，这对应了Helmholtz reciprocity原理），
    ray_position_image = [input_bubbletemp_silhouette(1,1),input_bubbletemp_silhouette(1,2),0];%c点 世界坐标系相机1图像点 发射光线
    d_maxValue_imagetobubble = abs(input_bubbletemp_silhouette(1,3));
    ray_direction_image = [0 0 -1];
    ray_image_world = [ray_position_image, ray_direction_image, d_maxValue_imagetobubble];   % c点 i1向量  随机采样，追踪到LED；i1在光锥角度omega_camera内
                
    %bubble_faceinfo_world = [input_bubbletemp_silhouette,input_bubbletemp_face];  %根据前面(3)的初步推断。  
    
    temp = dot(ray_image_world(1,4:6),input_bubbletemp_face);
    o1_temp_ray = ray_image_world(1,4:6)/rho_1_2-(temp+sqrt(abs(temp^2-(1-rho_1_2^2))))*...
        (input_bubbletemp_face/rho_1_2);          %Snell定律计算o1向量：
    % c+(d1i1+Δo1+d2o2)=q，->  d1<i1 x o2,o1> =<(q-c) x o2,o1>

    %推论1不共面定理，i1，o2和q-c共面的情况下，n1无法约束d1，
    raypoint_at_LED = input_bubbletemp_silhouette + dot(planeNormal_LED, planePoint_LED - input_bubbletemp_silhouette) / dot(planeNormal_LED, o1_temp_ray) * o1_temp_ray;
    %q点初步由o1与LED面交点求得
    tempcoplane = dot(cross(ray_direction_image,ray_direction_LED),o1_temp_ray);   %<i1×o2, o1>
    if tempcoplane ~= 0  %  if  <i1×o2, o1>≠0，此时推论1成立了，不共面定理
    %通过单根射线的追踪即可满足气泡前面的确定
        %raynumber_d1 = true; %则d1确定 end
         %根据d1<i1×o2, o1>=<(q-c)×o2, o1>，求得d1长度，和v1=c+d1i1点，
        d1_length_temp_ray = -dot(cross((raypoint_at_LED-ray_position_image),ray_direction_LED),o1_temp_ray)/tempcoplane; %求得d1长度即 d_imagetobubble
                            %这里的负号是因为是z轴的反方向的位置为视线方向。
        bubble_face_position = ray_position_image + d1_length_temp_ray*ray_direction_image; %v1=c+d1i1点，
    
        %代入o1=B(d1)[cosφ, sinφ]，B(d1)是q-v1和o2列空间的正交基；在欧几里德空间{R}}^{{3}}中，集合：{e1=(1,0,0), e2=(0,1,0), e3=(0,0,1)}组成一个标准正交基。
        %构建标准正交基B_d1
        v1 = (raypoint_at_LED-bubble_face_position)/norm((raypoint_at_LED-bubble_face_position));
        v2 = ray_direction_LED/norm(ray_direction_LED);%v3 = cross(v1,v2); v3 = v3/norm(v3);
        B_d1 = [v1' v2'];
        for degree_circle = 1:360
            phi = degree_circle*pi/180;
            face_direction = B_d1*[cos(phi);sin(phi)];  %这里求得一个待定的o1，需要进行判别
            %n⊥单位范数向量，与i1和n1共面，且正交于n1，基于snell定律μ1<i1, n ⊥ >=μ2<o1, n ⊥ >
            %μ1<i1-μ2/μ1·o1, n⊥ >=0，   i1-μ2/μ1·o1和n⊥垂直，且n1与n⊥垂直，所以i1-μ2/μ1·o1和n1平行
            n_parallel = ray_direction_image - rho_1_2*face_direction';
            n1 = input_bubbletemp_face;
            if dot(n1,n_parallel) <= 0.001 && dot(n1,n_parallel) >= -0.001 %判断n1与i1-μ2/μ1·o1是否平行
                %n1正确，d1求出，v1求出（这里会求出多个n1；气泡的正面信息）
                break;
            end
        end
        o1_tempupdate_ray = face_direction';
        n1_update = n_parallel;
        output_bubbletemp_face = ...
            [n1_update(1,1),n1_update(1,2),abs(n1_update(1,3))];%abs(n1_update(1,3)修正后的n1应该是朝向气泡外侧的        
        output_bubbletemp_silhouette = [input_bubbletemp_silhouette(1,1),input_bubbletemp_silhouette(1,2),d1_length_temp_ray];
    else
        output_bubbletemp_face = input_bubbletemp_face;
        output_bubbletemp_silhouette = input_bubbletemp_silhouette;
    end

end
end

%这个circle在q-v1和o2构成的正交基上，就是q-v1和o2所在平面上。
%{
aaap = zeros(360,3);
for n = 1:360
    phi = n*pi/180;
    face_direction = B_d1*[cos(phi);sin(phi)];
    aaap(n,1:3) = face_direction';
    plot3(face_direction(1,1),face_direction(2,1),face_direction(3,1),'ro')
    hold on
end
%}



