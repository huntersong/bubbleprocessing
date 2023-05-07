function [intersectionPoint] = intersectLinePlane(point,vector,planeNormal,planePoint)
% 定义直线和平面
%{
point = [1 2 3]; % 直线上的一个点
vector = [4 5 6]; % 直线的方向向量
planeNormal = [0 0 1]; % 平面的法向量
%}
% 计算平面上的一点
%planePoint = [0 0 0]; % 平面上的一个点
t = dot(planeNormal, planePoint - point) / dot(planeNormal, vector);
intersectionPoint = point + t * vector;
end


