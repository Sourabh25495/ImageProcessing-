clc;
clear all;
close all;
[filename, user_canceled] = imgetfile;
iout_hsv=rgb2hsv(im2double(imread(filename)));  
s=size(iout_hsv);
size1=s(1)*s(2);

%RMS contrast calc
total=sum(sum(iout_hsv(:,:,3)));
    imean=total/size1;
    sum1=0;
    
    
   for x=1:s(1)
        for y=1:s(2)
            b=(iout_hsv(x,y,3)-imean)^2;
            sum1=sum1+b;
            
            
           
            
        end
   end
    sum2=sum1/size1;
    cr=sqrt(sum2);
    cr
     
E = entropy(iout_hsv(:,:,3));
E
% Avg_E_local = mean(mean(entropyfilt(iout_hsv(:,:,3))))
%      
% Avg_sat=mean(mean(iout_hsv(:,:,2)))
