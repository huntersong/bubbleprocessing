# bubbleprocessing  
#ray tracing 程序分三个个部分，第一个部分是标定过程，获得相机和测量域在世界坐标系下的相对位置，本程序不涉及，直接采用正交的标定外参和内参；第二部分是预处理图像，对两个个相机的图像进行处理，获得光线追踪算法的输入初值；第三部分是基于光线追踪的三维重建  
使用方法：通过对当前文件夹里编辑image_camera1 = imread("figure1.bmp");选择本文件夹下的图片  
#ray tracing program is divided into three parts, the first part is the calibration process, to obtain the relative position of the camera and the measurement domain in the world coordinate system, this program is not involved, directly using orthogonal calibration external parameters and internal parameters; The second part is the preprocessed image, which processes the images of the two cameras to obtain the initial input values of the ray tracing algorithm. The third part is 3D reconstruction based on ray tracing  
Use: by editing image_camera1 = imread("figure1.bmp") in the current folder; Select the images in this folder  
![image](https://github.com/huntersong/bubbleprocessing/blob/main/bubbleimagegit/51bubble.png)
