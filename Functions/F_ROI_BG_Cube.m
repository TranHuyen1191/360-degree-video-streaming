%T(i) i=1:No_tile is selected versions for tile i
function [s_t_,T,sel_BR] = F_ROI_BG_Cube(angle,Fh, Fv, vp_W,vp_H, face_W, face_H, phi,theta,No_face,tile_hori_num,tile_ver_num,No_ver,Bw,BR,deltaBR,PSNR,MSE,deltaMSE,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Viewport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s_t_,s_,f_,m_,n_,m_c,n_c] = F_ExtractCubeTileCodOfVP(Fh, Fv, vp_W,vp_H, face_W, face_H, phi,theta,tile_hori_num,tile_ver_num,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H);

T = zeros(No_face,tile_hori_num*tile_ver_num);
Min_bw = 0; % BW min to transmit all face with lowest versions
sel_BR =0;
% if BW is enough --> transmit all face
for i = 1:No_face
	for k = 1:tile_hori_num*tile_ver_num 
		Min_bw = Min_bw + deltaBR(i,k,1);
	end
end

if Bw >= Min_bw 
	T = ones(No_face,tile_hori_num*tile_ver_num);
	sel_BR = Min_bw;
end

flag_stop = 0;
while flag_stop == 0
	flag_stop = 1; 
	for i = 1:No_face
		for k = 1:tile_hori_num*tile_ver_num 
            if s_t_(i,k)>0
                if T(i,k)+1 <= No_ver
                    temp_BR = sel_BR + deltaBR(i,k,T(i,k)+1);
                    if (temp_BR <= Bw)
                        T(i,k) = T(i,k)+1;
                        sel_BR = temp_BR;
                        flag_stop = 0; 
                    end
                end
            end
		end
	end
end





