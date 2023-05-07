
%作气泡前面的图bubble_all_position_face_world
bubble_front_position_face_world = zeros(length(bubble_all_position_face_world)/2,3);
bubble_front_face_direction = zeros(length(bubble_all_position_face_world)/2,3); 
for n=1:length(tempab)
    a=tempab(n,2)-tempab(n,1);
    b = (tempab(n,1)+1)/2;
    bubble_front_position_face_world(b:(b+a),1:3) = bubble_all_position_face_world(tempab(n,1):tempab(n,2),1:3);
    bubble_front_face_direction(b:(b+a),1:3) = bubble_all_face_direction(tempab(n,1):tempab(n,2),1:3);
end
figure(4)
ptCloud = pointCloud(bubble_front_position_face_world(:,1:3));   ptCloud;     %从输入点坐标创建一个点云对象。检查点云对象的属性
pcshow(ptCloud);
normals = pcnormals(ptCloud);   ptCloud = pointCloud(bubble_front_position_face_world(:,1:3),'Normal',normals); %在每个点上增加面法向：pcnormals评估这个点云上每个点的法向量
pcshow(ptCloud);
U_bubble_face_world = ptCloud.Location(:,1);    V_bubble_face_world = ptCloud.Location(:,2);...
    W_bubble_face_world = ptCloud.Location(:,3);  %点坐标
hold on;
quiver3(U_bubble_face_world,V_bubble_face_world,W_bubble_face_world,...
    bubble_front_face_direction(:,1),bubble_front_face_direction(:,2),bubble_front_face_direction(:,3),0.5);
hold off
view(0,180)

%作面的图bubble_silhouette
bubble_front_position_face_world = zeros(length(bubble_silhouette)/2,3);
bubble_front_face_direction = zeros(length(bubble_silhouette)/2,3); 
for n=1:length(tempab)
a=tempab(n,2)-tempab(n,1);
b = (tempab(n,1)+1)/2;
bubble_front_position_face_world(b:(b+a),1:3) = bubble_silhouette(tempab(n,1):tempab(n,2),1:3);
bubble_front_face_direction(b:(b+a),1:3) = bubble_all_face_direction(tempab(n,1):tempab(n,2),1:3);
end
figure(4)
ptCloud = pointCloud(bubble_front_position_face_world(:,1:3)); ptCloud; %从输入点坐标创建一个点云对象。检查点云对象的属性
pcshow(ptCloud);
normals = pcnormals(ptCloud); ptCloud = pointCloud(bubble_front_position_face_world(:,1:3),'Normal',normals); %在每个点上增加面法向：pcnormals评估这个点云上每个点的法向量
pcshow(ptCloud);
U_bubble_face_world = ptCloud.Location(:,1); V_bubble_face_world = ptCloud.Location(:,2);...
W_bubble_face_world = ptCloud.Location(:,3); %点坐标
hold on;
quiver3(U_bubble_face_world,V_bubble_face_world,W_bubble_face_world,...
bubble_front_face_direction(:,1),bubble_front_face_direction(:,2),bubble_front_face_direction(:,3),0.5);
hold off
view(0,180)

axis vis3d  %3维坐标系
for i = 1:36
    camorbit(10,0,'data',[0 1 0]) %%[0 0 1]表示按z轴旋转。36*10=360表示旋转一周
    drawnow %%即时显示旋转的结果
end
%https://blog.csdn.net/wayne6515/article/details/113980952
for i=1:36
    camorbit(10,0,'data',[0,1,0])%[0 0 1]表示按z轴旋转。36*10=360表示旋转一周
    M=getframe(gcf);
    nn=frame2im(M);
    [nn,cm]=rgb2ind(nn,256);
    if i==1
        imwrite(nn,cm,'figure3_out.gif','gif','LoopCount',inf,'DelayTime',0.1);%说明loopcount只是在i==1的时候才有用
    else
        imwrite(nn,cm,'figure3_out.gif','gif','WriteMode','append','DelayTime',0.1)%当i>=2的时候loopcount不起作用
    end
end




nn=1
for n=1:length(bubble_front_position_face_world)
    if bubble_front_position_face_world(n,3)>0
    bubble_front_position_face_world_new(nn,1:3) = bubble_front_position_face_world(n,1:3);
    nn = nn+1;
    end
end



