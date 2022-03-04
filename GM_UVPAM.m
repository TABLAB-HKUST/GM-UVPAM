%% GM scanning for UV-PAM

clear; clc;
folder = 'G:\LXF\UV-PAM\JoVE\'; 
filename = 'Example_1.dat';

%% 30 um scanning interval and 22 point
yn = 20000;     %% 1 for 0.3125 micronmeter.
zn = 128;        %% Smapling points on A-line signal
GMnum = 22;       %% Data point on 1/4 period of GM
xn = 168;         %% number of sub-iamges
cr=14;   %% calibration for starting of GM, may change from device from device
deltUp=2;  
deltD=2;  

%% GM scanning trace
Lx = 15; % half of the GM scanning interval
dy = 1*0.625/2; % step size in y axis
ddy = dy/GMnum;
c = yn*GMnum; %% total number of signal on one sub-image
t = 1:c;
yi = ddy*(t-1);
xi = Lx*sawtooth(2*pi*(yi/(4*dy)),0.5);
% figure, plot(yi,xi);

%% Extracting data from the measured data
fdaq = fopen(strcat(folder, filename));
for i = 1:xn
        daqall = fread(fdaq, zn*yn*GMnum, 'uint16');
        daqall = reshape(daqall, zn, yn*GMnum);
        pa1 = daqall-32767;   
        PA_MAP=max(pa1)-min(pa1); %% MAP

        %% calibration for the GM when it begains to scan
        PA_MAP(1:GMnum+cr)=[];  
        PA_MAP(end-2*GMnum+(GMnum+cr)+1:end)=[];

        PA_MAPall(i,:)=PA_MAP; %%
end
fclose(fdaq);

%% GM scanning trajectory
for i=1:xn
    x(i,:)=xi+Lx*(i-1); % xi is the sawtooth function, Lx is amplitude
    y(i,:)=yi;
	if mod(i,2)==0
        y(i,:)=fliplr(y(i,:));
	end
end

%% Reconstruct sub-images and stitching them together
Imgz=[];
for ii=1:xn
    z=PA_MAPall(ii,:);
    yy=y(ii,:);
    xx=x(ii,:);
    [Xi,Yi,Zi] = griddata(yy,xx,z,linspace(min(yy),max(yy), round((max(yy)-min(yy))/0.5))',linspace(min(xx),max(xx), (max(xx)-min(xx))/0.5+deltUp+deltD),'cubic');
    Zi=flipud(Zi);
    Zi(1:deltUp,:)=[];  %% delect the overlap point between adjacent sub-images
    Zi(end-deltD+1:end,:)=[];
    Imgz=[Zi;Imgz];
end
Imgz = medfilt2(Imgz,[3 3]);
figure;
imagesc(Imgz/max(max(Imgz)));
colormap(gray);
hold on;
axis equal;axis tight;

