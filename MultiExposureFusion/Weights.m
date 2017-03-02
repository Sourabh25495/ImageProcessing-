clc;
clear all;
close all;
n=input('enter number of images \n Eiffel Tower= 3 \n grandcanal=3 \n House=4  \n Igloo=6 \n Lincoln=4 \n Living Room=5 \n Lizard=9 \n Office=6 \n Studio=5 \n ');
ia=zeros(1,n);
iu=zeros();
io=zeros();
o=0;

for j=1:n
   [filename, user_canceled] = imgetfile;
i(:,:,:,j)=rgb2hsv(im2double(imread(filename)));  
end
% for j=1:n
% figure;imshow(i(:,:,:,j));
% end
s=size(i(:,:,:,1));
size1=s(1)*s(2);
for j=1:n
    total=sum(sum(i(:,:,3,j)));
    m=total/size1;
    if m>0.5
        ia(1,j)=1;
        o=o+1; %o gives number of underexposed images
    end
    u=n-o;%u gives number of underexposed images 
end

% calculation of weigts for underexposed
u2=2*u;
% u2=5;
sp=1/u2;%sp is the space between 2 x coordinates
x1(1,1)=sp/2;%x1 is the x coordinate according to number of images and x1(1,1) is the first cordinate
wu(1,1)=x1(1,1)*2;%w is the value of weights
for j=2:u2
x1(1,j)=x1(1,j-1)+sp;
x=x1(1,j);%x is temporary variable to compute values
if (x<=0.5)
    y=2*x;
elseif (0.5<x<=1)
    y=2*(1-x);
    
end
   wu(1,j)=y; 
end 
wusum=sum(sum(wu()));
wun=wu/wusum;%wun is the normalized weights for underexposed



% calculation of weigts for overexposed
o2=2*o;
% u2=5;
sp=1/o2;%sp is the space between 2 x coordinates
x1(1,1)=sp/2;%x1 is the x coordinate according to number of images and x1(1,1) is the first cordinate
wo(1,1)=x1(1,1)*2;%w is the value of weights
for j=2:o2
x1(1,j)=x1(1,j-1)+sp;
x=x1(1,j);%x is temporary variable to compute values
if (x<=0.5)
    y=2*x;
elseif (0.5<x<=1)
    y=2*(1-x);
    
end
   wo(1,j)=y; 
end 
wosum=sum(sum(wo()));
won=wo/wosum;%won is the normalized weights for over exposed
for j=(u+1):u2
    wun(1,j)=0;
end 
for j=1:o
    won(1,j)=0;
end

% for j=1:n
%     if (j<)
    w=[wun(1:u),won(o+1:o2)];
    fprintf('there are total %d images %d are under exposed %d are overexposed \n and the weights are as follows \n',n,u,o);
    w=w*2
for j=1:n
iw(:,:,:,j)=i(:,:,:,j)*w(1,j);%iw is the weighted image
end
for j=1:u
iu=iu+iw(:,:,:,j);
end
for j=(u+1):n
io=io+iw(:,:,:,j);
end
% i1=hsv2rgb(io);
% i2=hsv2rgb(iu);
% [h1 s1 v1] = rgb2hsv(i1);
% [h2 s2 v2] = rgb2hsv(i2);
v1= (io(:,:,3));
v2= (iu(:,:,3));
 m=rgb2hsv(io);
N = rgb2hsv(iu);
s=size(v1);
th1=0.95;
th2=0.05;
% fim=[];
for x=1:s(1)
    for y=1:s(2)
        ov=v1(x,y);
        un=v2(x,y);
        d=ov-un;
        if(un<th1)%black staurated
            f(x,y)=ov-d/4;
        elseif(ov>th2)
            f(x,y)=un+d/4;%white saturated 
        else
            f(x,y)=(((ov-d/4)*(1-ov))+((un+d/4)*un))/(1-d);
        end
%         hs(x,y)=h1(x,y)+s1(x,y)+f(x,y)+h2(x,y)+s2(x,y);   
% fim=[h1 s1 f(x,y)];
    end
end


% op=hsv2rgb(f);

% figure;imshow(v2),title('undersaturated image');
% figure;imshow(v1),title('oversaturated image');
% figure;imshow(f),title('output');
ihu=zeros();
isu=zeros();
iho=zeros();
iso=zeros();
for j=1:u
    ihu=i(:,:,1,j)+ihu;
%      isu=i(:,:,2,j)+isu;
end

ihau=ihu/u;%ihau is the hue averqage of all the underexposed images 


%   isau=isu/u;%isau is the saturation averqage of all the underexposed images 



for j=1+u:n
    iho=i(:,:,1,j)+iho;
%      iso=i(:,:,2,j)+iso;
end

ihao=iho/o;%ihau is the hue averqage of all the underexposed images 


%   isao=iso/o;%isau is the saturation averqage of all the underexposed images

%  % %% max satutration
for nn = 1 : u
    isao_all(:,:,nn) = i(:,:,2,nn);
    ihao_all(:,:,nn) = i(:,:,1,nn);
end
n_rem = n - u;
for nn = 1 : n_rem
    isau_all(:,:,nn) = i(:,:,2,nn+u);
    ihau_all(:,:,nn) = i(:,:,1,nn+u);
end
isao = mean(isao_all,3);
isau = mean(isau_all,3);
ihao = mean(ihao_all,3);
ihau = mean(ihau_all,3);

% isao = max(isao_all,[],3);
% isau = max(isau_all,[],3);
% ihao = max(ihao_all,[],3);
% ihau = max(ihau_all,[],3);

%  isao=zeros();
%  isau=zeros();
% for x=1:s(1)
% for y=1:s(2)
% isau=i(x,y,2,1);
% for j=2:u
%     if(i(x,y,2,j)>i(x,y,2,j-1))
%    isau(x,y,2,:)=i(x,y,2,j);   
%     end
% end
% end
% end
% 
% for x=1:s(1)
% for y=1:s(2)
% isao=i(x,y,2,u+1);
% for j=u+1:n
%     if(i(x,y,2,j)>i(x,y,2,j-1))
%    isao(x,y,2)=i(x,y,2,j);   
%     end
% end
% end
% end
% % isau=isu;
% % isao=iso;
% % 

%%
   
% %% kartalov color for 2 images
% 



h1=ihao;
h2=ihau;
s1=isao;
s2=isau;

for x=1:s(1)
    for y=1:s(2)
        sov(x,y)=s1(x,y);
        hov(x,y)=h1(x,y);
        sun(x,y)=s2(x,y);
        hun(x,y)=h2(x,y);
        vsfov(x,y)=s1(x,y)*(1-f(x,y));
         vsfun(x,y)=s2(x,y)*f(x,y);
         if((vsfov(x,y))<(vsfun(x,y)))
       sfu(x,y)=vsfun(x,y);
       hfu(x,y)=hun(x,y);
         else
              sfu(x,y)=sov(x,y);
       hfu(x,y)=hov(x,y);
         end
    
 end
end
HSV_fused(:,:,1)=hfu;
HSV_fused(:,:,2)=sfu;
HSV_fused(:,:,3)=f;
op=hsv2rgb(HSV_fused);
%  figure;imshow(i1,[]),title('image 1');
% %  figure;imshow(i2,[]),title('image 2');
figure;imshow(op,[]);title('output');



