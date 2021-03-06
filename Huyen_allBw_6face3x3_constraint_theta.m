% viewport to CM converstion
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
face_W = 960;
face_H = 960;
phi =0; % Longitude  of viewport
theta = 0; %altitude of viewport
No_face = 6;
tile_hori_num = 3;
tile_ver_num = 3;
No_tile = tile_hori_num * tile_ver_num;
No_frame = 300;

for j=1:tile_ver_num
    for i=1:tile_hori_num
        tile_id = (j-1) * tile_hori_num + i;
        LB_tile_W(tile_id) = face_W/tile_hori_num * (i-1);
        HB_tile_W(tile_id) = face_W/tile_hori_num * i;
        LB_tile_H(tile_id) = face_H/tile_ver_num * (j-1);
        HB_tile_H(tile_id) = face_H/tile_ver_num * j;
        
    end
end

               
No_ver = 9;
QP_ar = [50 48 44 40 36 32 28 24 20];
% tile version bitrates
BR = ones(No_face,No_tile,No_ver+1);

% tile version MSE
MSE = ones(No_face,No_tile,No_ver+1);
PSNR = zeros(No_face,No_tile,No_ver+1);
for f = 1:No_face		
	for t = 1: No_tile
       BR(f,t,:) = BR(f,t) .* zeros(1,No_ver+1);
       PSNR(f,t,:) = PSNR(f,t) .* zeros(1,No_ver+1);
       MSE(f,t,:) = MSE(f,t) .*  [65025 zeros(1,No_ver)];
    end
end


fname = sprintf('BR_PSNR_6f%.0fx%.0f_low_delay_%.0fFr.txt',tile_hori_num,tile_ver_num,No_frame);
fileID = fopen(fname,'r');
formatSpec = '%d\t%d'; % face(1) Tile(2)
for i =  1: No_ver 
    % Bitrate_encode(3) PSNR(4) MSE(5)                    
    formatSpec = strcat(formatSpec,'\t%f\t%f\t%f');
end
No_col = 2+3*No_ver;
sizeA = [No_col Inf];
A = fscanf(fileID,formatSpec,sizeA);

for i = 1:No_face
    for k = 1:No_tile
        index = (i-1)*No_tile+k;
        for j = 1:No_ver
           % PSNR (face,tile, ver+1)
           BR(i,k,j+1)= A(3+(j-1)*3,index);
           PSNR(i,k,j+1)= A(4+(j-1)*3,index);
           MSE(i,k,j+1)= A(5+(j-1)*3,index);
        end
    end
end



% calculate deltaMSE, deltaBR
deltaBR = zeros(No_face,No_tile,No_ver);
deltaMSE = zeros(No_face,No_tile,No_ver);
for i=1:No_ver
    deltaBR(:,:,i) = BR(:,:,i+1) - BR(:,:,i);
    deltaMSE(:,:,i) = MSE(:,:,i) - MSE(:,:,i+1);
end
   
fname = sprintf('Log%sH_Al_log_%.0fface%.0fx%.0f_constraint_vbr_theta.txt',bracket,No_face,tile_hori_num,tile_ver_num);
fout = fopen(fname,'w');

fprintf(fout,'PSNR_est\tT_e\tangle_theta\t');
for i = 1: No_face
    for k = 1: No_tile
      fprintf(fout,'Vf%dt%d\t',i,k);
    end
end

fprintf(fout,'curT_e\t');
for i = 1: No_face
    for k = 1: No_tile
      fprintf(fout,'Sf%dt%d\t',i,k);
    end
end
fprintf(fout,'\n');



BW_arr = [2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000];
for i = 1:1:1 %% Original 1:1:11
    T_e = BW_arr(i);
    for angle_phi = 0:30:0 
        if(angle_phi <=180)
            phi = angle_phi/180*pi;
        else
            phi = -(360-angle_phi)/180*pi;
        end
		
		for angle_theta = -165:15:-165  %% Original -165:15:180
			theta = angle_theta/180*pi;
			if(angle_theta >=0)
				angle_theta_name = angle_theta;
			else
				angle_theta_name = 360+angle_theta;
            end		
            %tic
			[TD_sel,v_sel,TB_sel,N_ft] = F_H_Adapt_constraint(Fh, Fv, vp_W,vp_H, face_W, ...
                face_H, phi,theta,No_face,tile_hori_num,tile_ver_num,No_ver,T_e,BR,...
                MSE,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H);
            %toc
            PSNR_sel = 10*log10(255*255/TD_sel);
			fprintf(fout,'%.0f\t%.0f\t%.0f\t',PSNR_sel,T_e,angle_theta_name);
			for i = 1: No_face
				for k = 1: No_tile
					fprintf(fout,'%.0f\t',v_sel(i,k)-1);
				end
			end

			fprintf(fout,'%.0f\t',TB_sel);
			for i = 1: No_face
				for k = 1: No_tile
					fprintf(fout,'%.0f\t',N_ft(i,k));
				end
			end
			fprintf(fout,'\n');
		end
    end
end
   
fclose all;