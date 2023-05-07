function [point1_img_camera1,point2_img_camera1] = slice_position(image) 
stats = regionprops('table', image, 'BoundingBox','Extrema');
boundbox = stats.BoundingBox;
%a = round(boundbox(1,1));
b = round(boundbox(1,2));
%width = round(boundbox(1,3))-1;
height = round(boundbox(1,4))-1;
point1_img_camera1 = zeros(height,3);
point2_img_camera1 = zeros(height,3);
for m=1:height
    mm = b+m;
    temp = regionprops('table', image(mm,:), 'BoundingBox','Extrema');
    aa = min(temp.BoundingBox(:,1));  aaa = sum(round(temp.BoundingBox(:,3)));
    point1_img_camera1(m,1) = round(aa);     point1_img_camera1(m,2) = mm;
    point2_img_camera1(m,1) = round(aa+aaa);     point2_img_camera1(m,2) = mm;
end
end
