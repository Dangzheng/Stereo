clear;clc;close all;

tex = imread('pic.jpg');
tex = double(rgb2gray(tex));
tex = tex/max(max(tex));% ��һ������
[x,y] = size(tex);
[x,y] = meshgrid(1:y,1:x);

%synthesize image size
width = 480;
height = 320;
depth = 225;
x_init = 2700;
y_init = 800;
% ����ͼ��Ĵ�С���ó�Ŀ��ͼ���2����������ͶӰ����ͼ���ʱ��
% ��ͶӰ����ͼ�����겻����������ô��������ͶӰ��ȥ�ĵ�������λ���ϲ�ֵ��
% ����������ͼ���ܹ���֤��ͶӰ��ͼ���ܹ��㹻�󸲸�right view��
[x_tex,y_tex] = meshgrid(x_init + 1:1:x_init + width*2,y_init + 1:1:y_init + height*2);
img_1 = interp2(x,y,tex,x_tex,y_tex);

%% ������ͼ������ϵ�µ�����
% ��ͼƬ��Ϊ������Ϣ��ʹ�ã���Ϊһ���������ƽ�档
[x1,y1] = meshgrid(1:1:width,1:1:height);
% �˴��Ƚ�Left view��ͼ��������
[x,y] = meshgrid(1:1:width*2,1:1:height*2);
img_l = interp2(x,y,img_1,x1,y1);

% �������������ͼ�������
x1 = reshape(x1.',1,[])';
y1 = reshape(y1.',1,[])';
%% �ڲξ������ξ���
%�ڲξ���
f = 325; % ����
xc = width/2;
yc = height/2;
    
intr_r=[f 0.0 xc;
        0.0 f yc;
        0.0 0.0 1.0];
intr_l=[f 0.0 xc;
        0.0 f yc;
        0.0 0.0 1.0];
%��x��Ϊ��ת��alpha�Ƕ�
% alpha = 20;
% R_r = [1,0,0,0;
%        0,cos(alpha*pi/180),sin(alpha*pi/180),0;
%        0,-sin(alpha*pi/180),cos(alpha*pi/180),0;
%        0,0,0,1];
%��y��Ϊ��ת��beta�Ƕ�
beta = 0;
R_r = [cos(beta*pi/180),0,-sin(beta*pi/180),0;
       0,1,0,0;
       sin(beta*pi/180),0,cos(beta*pi/180),0;
       0,0,0,1];
   
%Ϊ�˱�֤˫Ŀ֮��ı任��ϵ��������д��������������������ת��ת������
%Ȼ�����������ϵ���������ת������ֱ����R_l = R_related*R_r�õ�
%����������������֮��ֻ��ƽ��û���κε���ת
baseline = 2;  %��Ϊx��������Ϊ�ң����Դ˴���ֵӦ��Ϊ��ֵ��
R_related = [1,0,0,-baseline;
             0,1,0,0;
             0,0,1,0;
             0,0,0,1];
R_l = R_related*R_r;
% %% left view��right view��ͶӰ�����ڲ��ԵĴ��롣
% left_img = [x1';y1';ones(size(x1'))];% ת���������ʽ
% % ���������ϵ������
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
%% ����ͼƬ ����right view
[x_t,y_t] = meshgrid(1:1:width*2,1:1:height*2);
x_t = reshape(x_t.',1,[])';
y_t = reshape(y_t.',1,[])';
img_l1 = [x_t';y_t';ones(size(x_t'))];% ת���������ʽ
% ���������ϵ��
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

disp(['����������ɣ�'])
%% ����һ��ƽ���������
% ֻ��Ҫ�������һ��ƽ��Ȼ�󽻵��������õı���ͼƬ�Ͼ�ok�ˡ�
% �������������Ǿ����ȡһ��ͼƬ���ⲿ�ֵ�����ͼƬ��tri����ʾ

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
% right view ��β�ͬ��֮ǰ����ֻ��Ҫ����ͼ��С��ͼ����㹻�ˣ������ƶ�������������ģ����ǲ�care
x_tri = reshape(x_tri.',[],1)';
y_tri = reshape(y_tri.',[],1)';% �����ų���
left_img_tri = [x_tri;y_tri;ones(size(x_tri))];
left_cam_tri = intr_l^(-1)*left_img_tri;
[depth_tri,] = meshgrid(100:-0.1042:50,1:320);
depth_tri = reshape(depth_tri.',[],1)'; % �����ų���
left_cam_tri(4,:) = 1./depth_tri;
right_cam_tri = R_related*left_cam_tri;
right_img_tri = intr_r*right_cam_tri(1:3,:);
right_img_tri = right_img_tri(1:2,:)./repmat(right_img_tri(3,:),2,1);
img_pro = reshape(img_left_tri,1,[]);% ����չ���ų���
xq = reshape(right_img_tri(1,:),fliplr([height,width]))';xq = reshape(xq,1,[]);
yq = reshape(right_img_tri(2,:),fliplr([height,width]))';yq = reshape(yq,1,[]);
[x1,y1] = meshgrid(1:1:width,1:1:height);
img_right_tri = griddata(xq,yq,img_pro,x1,y1);
img_right_tri(isnan(img_right_tri)) = 0;
img_right_tri = reshape(img_right_tri,[height,width]);

%% �ϳ�
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
%% ����ground truth
% ������Ϊ����һ��б���һ��ƽ��ϳɵ�һ������ͼ������groundtruthҲҪ����������
% �����groundtruth����left view��ͼ��Ϊbase�����ġ�
depth_tri = reshape(depth_tri,fliplr([height,width]))'; % �����ų���
depth_tri(img_left_tri <= 0) = 0;
groundtruth = depth_tri;
groundtruth(depth_tri == 0) = depth; 
groundtruth = uint8(groundtruth);
figure(2);
imagesc(uint8(groundtruth)); f2 = getframe(gcf);colorbar;
imshow(f2.data);