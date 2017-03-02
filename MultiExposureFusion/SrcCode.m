
function varargout = sem8_gui(varargin)
% SEM8_GUI M-file for sem8_gui.fig
%      SEM8_GUI, by itself, creates a new SEM8_GUI or raises the existing
%      singleton*.
%
%      H = SEM8_GUI returns the handle to a new SEM8_GUI or the handle to
%      the existing singleton*.
%
%      SEM8_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEM8_GUI.M with the given input arguments.
%
%      SEM8_GUI('Property','Value',...) creates a new SEM8_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sem8_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sem8_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sem8_gui

% Last Modified by GUIDE v2.5 20-Jan-2016 16:45:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sem8_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sem8_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sem8_gui is made visible.
function sem8_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sem8_gui (see VARARGIN)

% Choose default command line output for sem8_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sem8_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sem8_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
choice=get(handles.popupmenu1,'value');
switch (choice)
 case 1
folder = 'House'; n = 4;


 case 2
folder = 'Eiffel_Tower'; n = 3;

 case 3
folder = 'Grandcanal'; n = 3;

%  case 4
% folder = 'Igloo'; n = 6;

 case 4
folder = 'Living_Room'; n = 5;

 case 5
folder = 'Studio'; n = 5;

 case 6
folder = 'Lizard'; n = 9;

 case 7
folder = 'Office'; n = 6;

 case 8
folder = 'Lincoln'; n = 4;
 
end
disp(folder);
disp(choice);

m=zeros(1,n);
o=0;

for j=1:n
    
    a = ['D:/cnds final/Exposure Fusion/' folder '/' num2str(j) '.jpg'];
    irgb(:,:,:,j)=(double(imread(a))/255);
    i(:,:,j)= rgb2gray(irgb(:,:,:,j));
    s=size(i(:,:,1));
end
save sem8_var.mat;


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load sem8_var.mat;
size1=s(1)*s(2);
for j=1:n
    total=sum(sum(i(:,:,j)));
    me=total/size1;
    m(1,j)=me;
    if me>0.5
        
        o=o+1; %o gives number of underexposed images
    end
    u=n-o;%u gives number of underexposed images
end

%%  simple avraging
isum=zeros();
irgb2=zeros();
irgb2=irgb;
for j=1:n
   
   isum=isum+irgb2(:,:,:,j);
end
isum=isum/n;
imshow(isum,[]),title('sa output');
save sem8_var.mat
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sem8_var.mat;
wsum=0;
sp=1/n;
x1(1,1)=sp/2;
w(1,1)=2*x1(1,1);
j=2;

while (j<=n)
    x1(1,j)=x1(1,j-1)+sp;
    x=x1(1,j);
    if(x<0.5)
        y=2*x;
    else
        y=2*(1-x);
    end
    w(1,j)=y;
     j=j+1;
end
for j=1:n
   wsum=wsum+w(1,j);
end
irgb1=zeros();
irgb1=irgb;
isum=zeros();
for j=1:n
    irgb1(:,:,:,j)=irgb1(:,:,:,j)*w(1,j);
   isum=isum+irgb1(:,:,:,j);
end

isum=isum/wsum;
imshow(isum,[]);title('wt output');
save sem8_var.mat;

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sem8_var.mat


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sem8_var.mat
switch (choice)
 case 1
% folder = 'House'; n = 4;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\house_tomm.jpg');title('tm output');


 case 2
folder = 'Eiffel_Tower'; n = 3;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\eiffel_tower tom.jpg');title('tm output');

 case 3
folder = 'Grandcanal'; n = 3;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\Grandcanal_tom.jpg');title('tm output');

%  case 4
% folder = 'Igloo'; n = 6;

 case 4
folder = 'Living_Room'; n = 5;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\Living_Room_tom.jpg');title('tm output');

 case 5
folder = 'Studio'; n = 5;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\Studio_tom.jpg');title('tm output');

 case 6
folder = 'Lizard'; n = 9;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\Lizard_tom.jpg');title('tm output');

 case 7
folder = 'Office'; n = 6;
imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\office_tom.jpg');title('tm output');

 case 8
folder = 'Lincoln'; n = 4;
 imshow('C:\Users\LNOVO\Desktop\tom_marten_output_images\lincoln_tom.jpg');title('tm output');

end 
save sem8_var.mat


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sem8_var.mat
for j=1:n
    
    a = ['D:/cnds final/Exposure Fusion/' folder '/' num2str(j) '.jpg'];
    irgb(:,:,:,j)=(double(imread(a))/255);
    i(:,:,j)= rgb2gray(irgb(:,:,:,j));
    s=size(i(:,:,1));
end
s=size(i(:,:,1));
size1=s(1)*s(2);
for j=1:n
    total=sum(sum(i(:,:,j)));
    me=total/size1;
    m(1,j)=me;
    if me>0.5
        
        o=o+1; %o gives number of underexposed images
    end
    u=n-o;%u gives number of underexposed images
end


%% calculatn of wts

sp=1/n;%sp is the space between 2 x coordinates
x1(1,1)=sp/2;%x1 is the x coordinate according to number of images and x1(1,1) is the first cordinate
w(1,1)=x1(1,1)*2;%w is the value of weights
for j=2:n
    x1(1,j)=x1(1,j-1)+sp;
    x=x1(1,j);%x is temporary variable to compute values
    if (x<=0.5)
        y=2*x;
    elseif (0.5<x<=1)
        y=2*(1-x);
        
    end
    w(1,j)=y;
end
wsum=sum(sum(w()));
wn=w/wsum;%wn is the normalized weights




%%
f=zeros();
for k =(u+1):n
    disp(k);
%     [h1 s1 v1] = rgb2hsv(i(:,:,:,k));
v1=i(:,:,k);
    % [h2 s2 v2] = rgb2hsv(iu);
    % v1= (io(:,:,3));
    % v2= (iu(:,:,3));
    % m=rgb2hsv(io);
    % N = rgb2hsv(iu);
    s=size(v1);
    
    % fim=[];
    for x=1:s(1)
        for y=1:s(2)
            ov=v1(x,y);
            if (m(1,k)<=ov<0.95)
                f(x,y,k)=ov*wn(1,k);
            elseif (0.05<ov< m(1,k))
                f(x,y,k)=ov;
            end
            
        end
    end
    
%     figure;imshow(i(:,:,k),[]);
%     figure;imshow(f(:,:,k),[]);
end


% f=zeros();
for k =1:u
    disp(k);
    v1=i(:,:,k);
%     
    s=size(v1);
    
   
    for x=1:s(1)
        for y=1:s(2)
            un=v1(x,y);
            if (m(1,k)>=un>0.05)
                f(x,y,k)=un*wn(1,k);
            elseif (0.95>un>m(1,k))
                f(x,y,k)=un;
            end
            
        end
    end
    
%     figure;imshow(i(:,:,k),[]);
%     figure;imshow(f(:,:,k),[]);
end

f = f + 1e-12;
fN = f ./ repmat(sum(f,3),[1 1 n]);
for k = 1 : n
%     figure; imshow(fN(:,:,k)); title('Norm weight');
end

%% Lets Call Tom Mertens algo

% create empty pyramid
pyr = gaussian_pyramid(zeros(s(1),s(2),3));
nlev = length(pyr);

% multiresolution blending
for i = 1:n
    % construct pyramid from each input image
	pyrW = gaussian_pyramid(fN(:,:,i));
	pyrI = laplacian_pyramid(irgb(:,:,:,i));
    
    % blend
    for l = 1:nlev
        w = repmat(pyrW{l},[1 1 3]);
        pyr{l} = pyr{l} + w.*pyrI{l};
    end
end

% reconstruct
R = reconstruct_laplacian_pyramid(pyr);
% figure;
 imshow(R);title('without postprocessing outpt');


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sem8_var.mat
for j=1:n
    
    a = ['D:/cnds final/Exposure Fusion/' folder '/' num2str(j) '.jpg'];
    irgb(:,:,:,j)=(double(imread(a))/255);
    i(:,:,j)= rgb2gray(irgb(:,:,:,j));
    s=size(i(:,:,1));
end
s=size(i(:,:,1));
size1=s(1)*s(2);
for j=1:n
    total=sum(sum(i(:,:,j)));
    me=total/size1;
    m(1,j)=me;
    if me>0.5
        
        o=o+1; %o gives number of underexposed images
    end
    u=n-o;%u gives number of underexposed images
end


%% calculatn of wts

sp=1/n;%sp is the space between 2 x coordinates
x1(1,1)=sp/2;%x1 is the x coordinate according to number of images and x1(1,1) is the first cordinate
w(1,1)=x1(1,1)*2;%w is the value of weights
for j=2:n
    x1(1,j)=x1(1,j-1)+sp;
    x=x1(1,j);%x is temporary variable to compute values
    if (x<=0.5)
        y=2*x;
    elseif (0.5<x<=1)
        y=2*(1-x);
        
    end
    w(1,j)=y;
end
wsum=sum(sum(w()));
wn=w/wsum;%wn is the normalized weights




%%
f=zeros();
for k =(u+1):n
    disp(k);
%     [h1 s1 v1] = rgb2hsv(i(:,:,:,k));
v1=i(:,:,k);
    % [h2 s2 v2] = rgb2hsv(iu);
    % v1= (io(:,:,3));
    % v2= (iu(:,:,3));
    % m=rgb2hsv(io);
    % N = rgb2hsv(iu);
    s=size(v1);
    
    % fim=[];
    for x=1:s(1)
        for y=1:s(2)
            ov=v1(x,y);
            if (m(1,k)<=ov<0.95)
                f(x,y,k)=ov*wn(1,k);
            elseif (0.05<ov< m(1,k))
                f(x,y,k)=ov;
            end
            
        end
    end
    
%     figure;imshow(i(:,:,k),[]);
%     figure;imshow(f(:,:,k),[]);
end


% f=zeros();
for k =1:u
    disp(k);
    v1=i(:,:,k);
%     
    s=size(v1);
    
   
    for x=1:s(1)
        for y=1:s(2)
            un=v1(x,y);
            if (m(1,k)>=un>0.05)
                f(x,y,k)=un*wn(1,k);
            elseif (0.95>un>m(1,k))
                f(x,y,k)=un;
            end
            
        end
    end
    
%     figure;imshow(i(:,:,k),[]);
%     figure;imshow(f(:,:,k),[]);
end

f = f + 1e-12;
fN = f ./ repmat(sum(f,3),[1 1 n]);
for k = 1 : n
%     figure; imshow(fN(:,:,k)); title('Norm weight');
end

%% Lets Call Tom Mertens algo

% create empty pyramid
pyr = gaussian_pyramid(zeros(s(1),s(2),3));
nlev = length(pyr);

% multiresolution blending
for i = 1:n
    % construct pyramid from each input image
	pyrW = gaussian_pyramid(fN(:,:,i));
	pyrI = laplacian_pyramid(irgb(:,:,:,i));
    
    % blend
    for l = 1:nlev
        w = repmat(pyrW{l},[1 1 3]);
        pyr{l} = pyr{l} + w.*pyrI{l};
    end
end

% reconstruct
R = reconstruct_laplacian_pyramid(pyr);
% figure;
% imshow(R);

% R1 = min(R(:));
% if R1 < 1
%     R = R - R1;
% end
% R =  R./max(R(:));
R (R<0)= 0;
R (R>1)= 1;
    
% figure; imshow(R);
iout_hsv=rgb2hsv(R);
iout_hsv(:,:,2)=iout_hsv(:,:,2)+0.1;
if iout_hsv(:,:,2)>1
    iout_hsv(:,:,2)=1;
end
% iout_hsv(:,:,3)=imadjust(iout_hsv(:,:,3),[0.1 0.9,[]]);
iout1=hsv2rgb(iout_hsv);
% figure;
imshow(iout1,[]);title('Proposed Algo pp output');

y = var(iout_hsv(:));
y
folder
