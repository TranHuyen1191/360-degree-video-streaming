% viewport to ERP converstion
clear;
clc
if ismac; bracket = '/'; else bracket = '\'; end;
addpath(genpath([pwd bracket 'Functions']))
addpath(genpath([pwd bracket 'Metadata']))
% Viewport
Fh = pi/2;
Fv = pi/2;
vp_W = 960;
vp_H = 960;
erp_W = 3840;
erp_H = 1920;
phi =0; % Longitude  of viewport
theta = 0; %altitude of viewport
tile_hori_num = 8;
tile_ver_num = 8;
No_tile = tile_hori_num * tile_ver_num;
No_frame = 300;

for j=1:tile_ver_num
    for i=1:tile_hori_num
        tile_id = (j-1) * tile_hori_num + i;
        LB_tile_W(tile_id) = erp_W/tile_hori_num * (i-1);
        HB_tile_W(tile_id) = erp_W/tile_hori_num * i;
        LB_tile_H(tile_id) = erp_H/tile_ver_num * (j-1);
        HB_tile_H(tile_id) = erp_H/tile_ver_num * j;
        
    end
end
                
No_ver = 9;
QP_ar = [50 48 44 40 36 32 28 24 20];
% tile version bitrates
BR = ones(No_tile,No_ver+1) .* zeros(1,No_ver+1);

% tile version MSE
MSE = ones(No_tile,No_ver+1) .* [65025 zeros(1,No_ver)];

fname = sprintf('BR_PSNR_tile_%.0fx%.0f_low_delay_%.0fFr.txt',tile_hori_num,tile_ver_num,No_frame);
fileID = fopen(fname,'r');
formatSpec = '%d'; % Tile(1)
for i =  1: No_ver 
    % Bitrate_encode(2) PSNR(3) MSE(4)                    
    formatSpec = strcat(formatSpec,'\t%f\t%f\t%f');
end
No_col = 1+3*No_ver;
sizeA = [No_col Inf];
A = fscanf(fileID,formatSpec,sizeA);

for i = 1:No_tile
    for j = 1:No_ver
        % BR (tile, ver+1)
        BR(i,j+1)= A(2+(j-1)*3,i);
        % MSE (tile, ver+1)
        MSE(i,j+1)= A(4+(j-1)*3,i);
    end
end

% calculate Utility and Cost
Cti = zeros(No_tile, No_ver); % Cost
Uti = zeros(No_tile, No_ver); % Utility
for i=1:No_ver
    Cti(:,i) = BR(:,i+1) - BR(:,i);
    Uti(:,i) = MSE(:,i) - MSE(:,i+1);
end
   
   
fname = sprintf('Log%sBellLab_log_tile%.0fx%.0f.txt',bracket,tile_hori_num,tile_ver_num);
fout = fopen(fname,'w');

fprintf(fout,'BW\tangle\t');
for i = 1: No_tile
  fprintf(fout,'V%d\t',i);
end

fprintf(fout,'curBW\t');
for i = 1: No_tile
  fprintf(fout,'P%d\t',i);
end
fprintf(fout,'\n');

BW_arr = [2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000];
for i = 1:1:1 %% Original 1:1:11
    Bw = BW_arr(i);
    for angle = 0:350:350 %% Original 0:15:350
        if(angle <=180)
            phi = angle/180*pi;
        else
            phi = -(360-angle)/180*pi;
        end
        [T,current_BW,m_,n_,P,B] = F_BellLab(Fh, Fv, vp_W,vp_H, erp_W, erp_H, phi,theta,No_tile,No_ver,Bw,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H,Uti,Cti);
        fprintf(fout,'%.0f\t%.0f\t',Bw,angle);
        for i = 1: No_tile
          fprintf(fout,'%.0f\t',T(i));
        end

        fprintf(fout,'%.0f\t',current_BW);
        for i = 1: No_tile
          fprintf(fout,'%.0f\t',P(i));
        end
        fprintf(fout,'\n');
    end
end
   
fclose all;
   