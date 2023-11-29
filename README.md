# bubbleprocessing  
该程序对单幅图像进行光线追踪，获得气泡的三维形态中的细节信息。
#ray tracing 程序分三个个部分：
第一个部分是标定过程，获得相机和测量域在世界坐标系下的相对位置以及可以获得气泡之间的相对空间信息，这里需要最基本的双目标定，本程序不涉及，直接采用正交布置的相机标定外参和内参；
第二部分是预处理图像，对两个个相机的图像进行处理，获得光线追踪算法的输入初值，这里主要获得单个气泡的二维图像；
第三部分是基于光线追踪的三维重建  
使用方法：运行raytrace_main.m程序
通过对当前文件夹里编辑image_camera1 = imread("figure1.bmp");选择本文件夹下的图片，或者其他图像，最理想的图像是单个气泡，气泡的背景颜色为白色，仅仅气泡区域存在灰度差异，便于直接进行光线追踪。  

This passage describes a process involving three main parts:

Calibration Process:
Obtaining the relative positions of the cameras and the measurement domain in the world coordinate system, as well as acquiring the relative spatial information between bubbles. This requires basic stereo calibration. The program uses the external and internal parameters of cameras arranged orthogonally without involving stereo calibration directly.

Image Preprocessing:
Processing the images from two cameras to obtain the initial values for the ray tracing algorithm. This step primarily involves obtaining the two-dimensional image of a single bubble.

Three-Dimensional Reconstruction Based on Ray Tracing:
This is the main part of the process, where three-dimensional reconstruction is performed using ray tracing.

Usage:
Run the "raytrace_main.m" program. Choose the image by editing the line "image_camera1 = imread("figure1.bmp");" in the current folder. The preferred image is one with a single bubble, where the background is white, and only the bubble region exhibits grayscale differences, facilitating direct ray tracing.
more url https://huntersong.com/2023/10/31/%e7%89%b9%e6%ae%8a%e7%9a%84%e6%b0%94%e6%b3%a1/  
![image](https://github.com/huntersong/bubbleprocessing/blob/main/bubbleimagegit/51bubble.png)
