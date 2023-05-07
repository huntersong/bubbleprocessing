function [Point_camera1and2] = slice_camera1and2_position(point1,point2,camera1,camera2) 
    % 两条直线的参数形式：a/bx + z = c/b 
    ab1 = tan(camera1); ab2 = -tan(camera2);
    cb1 = ab1*point1(1,1)+0*point1(1,3);
    cb2 = ab2*point2(1,1)+1*point2(1,3);
    line1 = [ab1, 0];
    line2 = [ab2, 1];
    Point = [line1;line2]\[cb1;cb2];
    Point_camera1and2 = [Point(1,1),point1(1,2),Point(2,1)];
end

