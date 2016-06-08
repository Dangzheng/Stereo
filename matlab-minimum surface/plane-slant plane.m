clear;clc;close all;

tex = imread('pic.jpg');
tex = double(rgb2gray(tex));
tex = tex/max(max(tex));% 归一化处理
[x,y] = size(tex);
[x,y] = meshgrid(1:y,1:x);

%synthesize image size
width = 480;
height = 320;
depth = 225;
x_init = 2700;
y_init = 800;
% 纹理图像的大小设置成目标图像的2倍，这样重投影到右图像的时候
% 重投影到右图的坐标不是整数，那么就利用重投影回去的点在整数位置上插值，
% 两倍的纹理图像能够保证重投影的图像能够足够大覆盖right view。
[x_tex,y_tex] = meshgrid(x_init + 1:1:x_init + width*2,y_init + 1:1:y_init + height*2);
img_1 = interp2(x,y,tex,x_tex,y_tex);

%% 创建左图像坐标系下的坐标
% 将图片作为纹理信息来使用，作为一个带纹理的平面。
[x1,y1] = meshgrid(1:1:width,1:1:height);
% 此处先将Left view的图像做出来
[x,y] = meshgrid(1:1:width*2,1:1:height*2);
img_l = interp2(x,y,img_1,x1,y1);

% 这里的坐标是左图像的坐标
x1 = reshape(x1.',1,[])';
y1 = reshape(y1.',1,[])';
%% 内参矩阵和外参矩阵
%内参矩阵
f = 325; % 焦距
xc = width/2;
yc = height/2;
    
intr_r=[f 0.0 xc;
        0.0 f yc;
        0.0 0.0 1.0];
intr_l=[f 0.0 xc;
        0.0 f yc;
        0.0 0.0 1.0];
%以x轴为轴转动alpha角度
% alpha = 20;
% R_r = [1,0,0,0;
%        0,cos(alpha*pi/180),sin(alpha*pi/180),0;
%        0,-sin(alpha*pi/180),cos(alpha*pi/180),0;
%        0,0,0,1];
%以y轴为轴转动beta角度
beta = 0;
R_r = [cos(beta*pi/180),0,-sin(beta*pi/180),0;
       0,1,0,0;
       sin(beta*pi/180),0,cos(beta*pi/180),0;
       0,0,0,1];
   
%为了保证双目之间的变换关系，我们先写出他们在右相机向左相机转换转换矩阵
%然后从世界坐标系到左相机的转换矩阵直接由R_l = R_related*R_r得到
%假设左相机到右相机之间只有平移没有任何的旋转
baseline = 2;  %因为x轴正方向为右，所以此处的值应该为负值；
R_related = [1,0,0,-baseline;
             0,1,0,0;
             0,0,1,0;
             0,0,0,1];
R_l = R_related*R_r;
% %% left view向right view的投影，用于测试的代码。
% left_img = [x1';y1';ones(size(x1'))];% 转换成齐次形式
% % 左相机坐标系下坐标
% left_cam = intr_l^(-1)*left_img;
% left_cam(4,:)= 1/depth * ones(size(x1'));
% right_cam = R_related*left_cam;
% right_img = intr_r*right_cam(1:3,:);
% right_img = right_img(1:2,:)./repmat(right_img(3,:),2,1);
% img_pro = reshape(img_l,1,[]);
% xq = reshape(right_img(1,:),fliplr([height,width]))';xq = reshape(xq,1,[]);
% yq = reshape(right_img(2,:),fliplr([height,width]))';yq = reshape(yq,1,[]);
% x1 = reshape(x1,fliplr([height,width]))';
% y1 = reshape(y1,fliplr([height,width]))';
% img_r = griddata(xq,yq,img_pro,x1,y1);
% img_r(isnan(img_r)) = 0;
% img_r = reshape(img_r,[height,width]);
% figure(2);
% imshow(img_r);
%% 纹理图片 制作right view
[x_t,y_t] = meshgrid(1:1:width*2,1:1:height*2);
x_t = reshape(x_t.',1,[])';
y_t = reshape(y_t.',1,[])';
img_l1 = [x_t';y_t';ones(size(x_t'))];% 转换成齐次形式
% 左相机坐标系下
left_cam = intr_l^(-1)*img_l1;
left_cam(4,:)= 1/depth * ones(size(x_t'));
right_cam = R_related*left_cam;
right_img = intr_r*right_cam(1:3,:);
right_img = right_img(1:2,:)./repmat(right_img(3,:),2,1);
img_pro = reshape(img_1,1,[]);
xq = reshape(right_img(1,:),fliplr([height*2,width*2]))';xq = reshape(xq,1,[]);
yq = reshape(right_img(2,:),fliplr([height*2,width*2]))';yq = reshape(yq,1,[]);
x1 = reshape(x1,fliplr([height,width]))';
y1 = reshape(y1,fliplr([height,width]))';
img_r = griddata(xq,yq,img_pro,x1,y1);
img_r(isnan(img_r)) = 0;
img_r = reshape(img_r,[height,width]);

disp(['背景制作完成！'])
%% 另外一个平面的制作：
% 只需要随便制作一个平面然后交叠在制作好的背景图片上就ok了。
% 首先先用上三角矩阵截取一个图片。这部分的纹理图片用tri来表示

tri = triu(ones(height*2,height*2),0);
x_init_tri = 2800;
y_init_tri = 1750;
[x,y] = size(tex);
[x,y] = meshgrid(1:y,1:x);
[x_tex_tri,y_tex_tri] = meshgrid(1+x_init_tri:1:height*2+x_init_tri,1+y_init_tri:1:height*2+y_init_tri);
img_2 = interp2(x,y,tex,x_tex_tri,y_tex_tri);
% imshow(img_2);
img_2(tri>0) = 0;
[x_tex_tri,y_tex_tri] = size(img_2);
[x_tex_tri,y_tex_tri] = meshgrid(1:y_tex_tri,1:x_tex_tri);
[x_tri,y_tri] = meshgrid(1:1:width,1:1:height);
% left view
img_left_tri = interp2(x_tex_tri,y_tex_tri,img_2,x_tri,y_tri+100);
% right view 这次不同于之前的是只需要用左图大小的图像就足够了，他爱移动到哪里是随意的，我们不care
x_tri = reshape(x_tri.',[],1)';
y_tri = reshape(y_tri.',[],1)';% 按行排成行
left_img_tri = [x_tri;y_tri;ones(size(x_tri))];
left_cam_tri = intr_l^(-1)*left_img_tri;
[depth_tri,] = meshgrid(100:-0.1042:50,1:320);
depth_tri = reshape(depth_tri.',[],1)'; % 按行排成行
left_cam_tri(4,:) = 1./depth_tri;
right_cam_tri = R_related*left_cam_tri;
right_img_tri = intr_r*right_cam_tri(1:3,:);
right_img_tri = right_img_tri(1:2,:)./repmat(right_img_tri(3,:),2,1);
img_pro = reshape(img_left_tri,1,[]);% 按列展开排成列
xq = reshape(right_img_tri(1,:),fliplr([height,width]))';xq = reshape(xq,1,[]);
yq = reshape(right_img_tri(2,:),fliplr([height,width]))';yq = reshape(yq,1,[]);
[x1,y1] = meshgrid(1:1:width,1:1:height);
img_right_tri = griddata(xq,yq,img_pro,x1,y1);
img_right_tri(isnan(img_right_tri)) = 0;
img_right_tri = reshape(img_right_tri,[height,width]);

%% 合成
img_l(img_left_tri > 0) = 0;img_l = img_l + img_left_tri;
img_r(img_right_tri > 0) = 0;img_r = img_r + img_right_tri;
figure(1);
scr = get(0,'ScreenSize');
set(gcf,'Position',[150,scr(4)/5,scr(3)-300,scr(4)/2]);
display = axes('Parent',gcf,... 
    'Units','pixels',...  
    'Position',[50 75 width height]); 
imshow(img_l);
hold on;
pause(0.01);
display = axes('Parent',gcf,... 
    'Units','pixels',...  
    'Position',[scr(3)- 300-(width+50),75,width,height]); 
imshow(img_r);
%% 制作ground truth
% 这里因为是用一个斜面和一个平面合成的一个虚拟图像，所以groundtruth也要分两步来做
% 这里的groundtruth是以left view的图像为base制作的。
depth_tri = reshape(depth_tri,fliplr([height,width]))'; % 按行排成行
depth_tri(img_left_tri <= 0) = 0;
groundtruth = depth_tri;
groundtruth(depth_tri == 0) = depth; 
groundtruth = uint8(groundtruth);
figure(2);
imagesc(uint8(groundtruth)); f2 = getframe(gcf);colorbar;
imshow(f2.data);