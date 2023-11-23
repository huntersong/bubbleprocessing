function showfigure(bubble_all_position_face_world)
figure(1)
ptCloud = pointCloud(bubble_all_position_face_world(:,1:3));   ptCloud;     %从输入点坐标创建一个点云对象。检查点云对象的属性
pcshow(ptCloud);
%{
U_bubble_face_world = ptCloud.Location(:,1);    V_bubble_face_world = ptCloud.Location(:,2);...
    W_bubble_face_world = ptCloud.Location(:,3);  %点坐标
U_normals_bubble_face_world = bubble_all_face_direction(:,1);V_normals_bubble_face_world = bubble_all_face_direction(:,2);...
    W_normals_bubble_face_world = bubble_all_face_direction(:,3);
%U_normals_bubble_face_world = bubble_all_position_face_world(:,4);    
%V_normals_bubble_face_world = bubble_all_position_face_world(:,5);...
%W_normals_bubble_face_world = bubble_all_position_face_world(:,6); %构建气泡表面的法向量
hold on;
quiver3(U_bubble_face_world,V_bubble_face_world,W_bubble_face_world,...
    U_normals_bubble_face_world,V_normals_bubble_face_world,W_normals_bubble_face_world,0.5);
hold off
%}
view(0,180)
end