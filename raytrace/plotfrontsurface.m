
%获得气泡前面bubble_all_position_face_world
for numb=1:5:26
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'data.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'datatempab.mat'];
    load(filename)
    bubble_front_position_face_world = zeros(length(bubble_all_position_face_world)/2,3);
    bubble_back_position_face_world = zeros(length(bubble_all_position_face_world)/2,3);
    bubble_front_face_direction = zeros(length(bubble_all_position_face_world)/2,3); 
    for n=1:length(tempab)
        a = tempab(n,2)-tempab(n,1);
        b = (tempab(n,1)+1)/2;
        bubble_front_position_face_world(b:(b+a-1),1:3) = bubble_all_position_face_world(tempab(n,1):(tempab(n,2)-1),1:3);
        %bubble_front_face_direction(b:(b+a),1:3) = bubble_all_face_direction(tempab(n,1):tempab(n,2),1:3);
        bubble_back_position_face_world(b:(b+a-1),1:3) = bubble_all_position_face_world(tempab(n,2):(tempab(n,2)+a-1),1:3);
    end
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'datafrontface.mat'];
    %bubble_all_position_face_world = [bubble_all_position_face_world(:,2),bubble_all_position_face_world(:,1),...
    %    bubble_all_position_face_world(:,3),bubble_all_position_face_world(:,4:6)];
    save(filename,'bubble_front_position_face_world');
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'databackface.mat'];
    save(filename,'bubble_back_position_face_world');

%绘制光滑的面
    countmax = max(tempab(:,2)-tempab(:,1));
    [unique_nums, ~, idx] = unique(bubble_front_position_face_world(:,2));
    count = histcounts(idx, 1:numel(unique_nums)+1);
    positions = accumarray(idx(:), (1:numel(bubble_front_position_face_world(:,2))).', [], @(x) {sort(x)});
    
    bubble_position_X = zeros(length(count),countmax);
    bubble_position_Y = zeros(length(count),countmax);
    bubble_position_Z = zeros(length(count),countmax);
    %气泡的正面
    for n =1:length(count)
        num = length(count)-n+1;
        counttemp_a = min(str2num(num2str(positions{num })));    counttemp_b = max(str2num(num2str(positions{num })));
        countslice = count(num );  
        
        xyz = bubble_front_position_face_world(counttemp_a:counttemp_b,1:3);  % 原始数组
        x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
        % 插值扩展
        x_interp = linspace(1, numel(x), countmax);  % 生成插值后的横坐标
        % 进行插值
        bubble_position_x = interp1(x, x_interp, 'cubic');  % 使用线性插值方法
        bubble_position_y = interp1(y, x_interp, 'cubic');  % 使用线性插值方法
        bubble_position_z = interp1(z, x_interp, 'cubic');  % 使用线性插值方法
    
        bubble_position_X(n,:) = bubble_position_x;
        bubble_position_Y(n,:) = bubble_position_y;
        bubble_position_Z(n,:) = bubble_position_z;

    end
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'datafrontface_surf.mat'];
    save(filename,'bubble_position_X','bubble_position_Y','bubble_position_Z');


    tempabback(:,1) = tempab(:,2);
    tempabback(:,2) = tempab(:,2)+tempab(:,2)-tempab(:,1)-1;
    countmax = max(tempabback(:,2)-tempabback(:,1));    
    [unique_nums, ~, idx] = unique(bubble_back_position_face_world(:,2));
    count = histcounts(idx, 1:numel(unique_nums)+1);
    positions = accumarray(idx(:), (1:numel(bubble_back_position_face_world(:,2))).', [], @(x) {sort(x)});
      
    bubble_position_Xback = zeros(length(count),countmax);
    bubble_position_Yback = zeros(length(count),countmax);
    bubble_position_Zback = zeros(length(count),countmax);
    %气泡的背面
    for n =1:length(count)
        num = length(count)-n+1;
        counttemp_a = min(str2num(num2str(positions{num })));    counttemp_b = max(str2num(num2str(positions{num })));
        countslice = count(num );  
        
        xyz = bubble_back_position_face_world(counttemp_a:counttemp_b,1:3);  % 原始数组
        x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
        % 插值扩展
        x_interp = linspace(1, numel(x), countmax);  % 生成插值后的横坐标
        % 进行插值
        bubble_position_xx = interp1(x, x_interp, 'cubic');  % 使用线性插值方法
        bubble_position_yy = interp1(y, x_interp, 'cubic');  % 使用线性插值方法
        bubble_position_zz = interp1(z, x_interp, 'cubic');  % 使用线性插值方法
    
        bubble_position_Xback(n,:) = bubble_position_xx;
        bubble_position_Yback(n,:) = bubble_position_yy;
        bubble_position_Zback(n,:) = bubble_position_zz;

    end
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'databackface_surf.mat'];
    save(filename,'bubble_position_Xback','bubble_position_Yback','bubble_position_Zback');
    clear
end

%保存每个图
for numb = 1:5:26
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'datafrontface_surf.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'databackface_surf.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'datafrontface.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'databackface.mat'];
    load(filename)
    
    s = surf(bubble_position_X,bubble_position_Y,bubble_position_Z);
    colorbar
    colormap(gray(10))
    %sc = surfc(bubble_position_X,bubble_position_Y,bubble_position_Z,'FaceColor','b');
    hold on
    scatter3(bubble_back_position_face_world(:,1),bubble_back_position_face_world(:,2), ...
        bubble_back_position_face_world(:,3),'k.')
    colorbar;colormap(gray(10))
    %ptCloud = pointCloud(bubble_front_position_face_world(:,1:3));   ptCloud;    
    %pcshow(ptCloud);
    %surf(bubble_position_Xback,bubble_position_Yback,bubble_position_Zback)
    
    grid off % 去掉坐标网格 
    shading interp 
    axis equal
    axis tight
    s.FaceAlpha = 0.9;
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    view(10,80)
    %direction = [1 0 0];
    %rotate(s,direction,15)
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
    set(gcf,'color','w');
    xlim([min(min(bubble_position_X)) max(max(bubble_position_X))]);
    ylim([min(min(bubble_position_Y)) max(max(bubble_position_Y))]);
    zlim([min(min(bubble_back_position_face_world(:,3))) max(max(bubble_position_Z))]);
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'bubble.fig'];
    savefig(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'bubble'];
    saveas(gcf,filename,'tiffn')
    clear
    close all
end

%合并的图片
bh = 0;
for numb = 6:5:26
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'datafrontface_surf.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'databackface_surf.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'datafrontface.mat'];
    load(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\',sprintf('%02d',numb),'databackface.mat'];
    load(filename)
    bh = bh + 1;
    bubble_height = 250*(bh-1);
    s = surf(bubble_position_X,bubble_position_Y+bubble_height,bubble_position_Z);
    colorbar
    colormap(gray(10))
    %sc = surfc(bubble_position_X,bubble_position_Y,bubble_position_Z,'FaceColor','b');
    hold on
    scatter3(bubble_back_position_face_world(:,1),bubble_back_position_face_world(:,2)+bubble_height, ...
        bubble_back_position_face_world(:,3),'k.')
    colorbar;colormap(gray(10))
    %ptCloud = pointCloud(bubble_front_position_face_world(:,1:3));   ptCloud;    
    %pcshow(ptCloud);
    %surf(bubble_position_Xback,bubble_position_Yback,bubble_position_Zback)
    hold on
    grid off % 去掉坐标网格 
    shading interp 
    axis equal
    axis tight
    s.FaceAlpha = 0.9;
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    view(10,80)

    %direction = [1 0 0];
    %rotate(s,direction,15)
    set(gcf,'color','w');
    zlim([min(min(bubble_back_position_face_world(:,3))) max(max(bubble_position_Z))]);
    %{
    xlim([min(min(bubble_position_X)) max(max(bubble_position_X))]);
    ylim([min(min(bubble_position_Y)) max(max(bubble_position_Y))]);
    zlim([min(min(bubble_position_Z)) max(max(bubble_position_Z))]);
    
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'bubble.fig'];
    savefig(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case1\',sprintf('%02d',numb),'bubble'];
    saveas(gcf,filename,'tiffn')
    clear
    close all
    %}
end

set(gcf, 'position', [100 100 500 900]);
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\','bubbleall.fig'];
    savefig(filename)
    filename = ['D:\办公电脑文件\reynoldsstress\Motion Detection\CEJ\case6\','bubbleall.fig'];
    saveas(gcf,filename,'tiffn')



%绘制旋转的GIF图
%{
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
%}

%{
[x,y,z] = meshgrid([0:0.5:250],[140:0.5:320],[-1020:0.5:-880]);

V = zeros(361,501,281);
for xx=1:361
    for yy=1:501
        for zz=1:281
            for n=1:42642
            deth = sqrt((x(xx,yy,zz)-bubble_position_X3D(n,1))^2+(y(xx,yy,zz)-bubble_position_Y3D(n,1))^2+ ...
                (z(xx,yy,zz)-bubble_position_Z3D(n,1))^2);
                if deth < 1
                    V(xx,yy,zz) = 1;
                else
                    V(xx,yy,zz) = 0;
                end
            end
            fprintf("Iteration %d.\n", n);
        end
    end
end
%}
%{
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

x = bubble_front_position_face_world(:,1);
y = bubble_front_position_face_world(:,2);
freq = bubble_front_position_face_world(:,3);
tri = delaunay(x,y);
h = trisurf(tri, x, y, freq);
view(3);
axis vis3d;
lighting phong;
shading interp;
zlim([-1100 -800]) 

%}




