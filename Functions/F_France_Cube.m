%T(i) i=1:No_tile is selected versions for tile i
function [N_ft,v_sel,TB_sel,MSE_over_sel] = F_France_Cube(Fh, Fv, vp_W,vp_H, ...
    face_W, face_H, phi,theta,No_face,tile_hori_num,tile_ver_num,...
    No_ver,T_e,BR,MSE,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H)
rb = 3.5;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Viewport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N_ft,N_f,f_,m_,n_,m_c,n_c] = F_ExtractCubeTileCodOfVP(Fh, Fv, vp_W,vp_H, face_W, face_H, phi,theta,tile_hori_num,tile_ver_num,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the selected version for each tile
v_sel = ones(No_face,tile_hori_num*tile_ver_num); % array of selected versions
v_temp = ones(No_face,tile_hori_num*tile_ver_num); % initialize 1 (no transmit)
v_sel_vi = 1; % selected versions of visible tiles
v_sel_in = 1; % selected versions of invisible tiles
TB_sel =0;
TB_temp =0;
MSE_over_temp = 65025;
MSE_over_sel = 65025;
Min_br = 0; % T_e min to transmit all face with lowest versions
No_tile_vi = 0; % Number of visible tiles
No_tile_in = 0; % Number of invisible tiles
Aver_BR_vi = 0; % Average bitrate of visible tiles
Aver_BR_in = 0; % Average bitrate of invisible tiles
Sum_BR_vi_temp = 0; 
Sum_BR_in_temp = 0; 
No_tile = tile_hori_num*tile_ver_num;
N_ft_max = 0;

%% Calculate Min_br
for i = 1:No_face
	for k = 1: No_tile
		Min_br = Min_br + BR(i,k,2);
        if(N_ft(i,k) >N_ft_max)
            No_tile_vi = No_tile_vi+1;
            F_vi(No_tile_vi) = i;
            T_vi(No_tile_vi) = k;
        else
            No_tile_in = No_tile_in +1;
            F_in(No_tile_in) = i;
            T_in(No_tile_in) = k;
        end
	end
end

if T_e >= Min_br 
	v_sel = ones(No_face,tile_hori_num*tile_ver_num)*2;
    v_sel_vi = 2;
    v_sel_in = 2;
	TB_sel = Min_br;
end

%% Find selected versions for all tiles
v_temp = v_sel;
for v_in = v_sel_in:No_ver+1
    for cnt_in = 1:No_tile_in
        v_temp(F_in(cnt_in),T_in(cnt_in)) = v_in;
    end
    for v_vi = v_in : No_ver+1
        for cnt_vi = 1:No_tile_vi
            v_temp(F_vi(cnt_vi),T_vi(cnt_vi)) = v_vi;
        end
        TB_temp  = 0;
        Sum_BR_vi_temp = 0;
        Sum_BR_in_temp = 0;
        MSE_over_temp = 0;
        % Calculate TB
        for i = 1:No_face
            for k = 1:tile_hori_num*tile_ver_num 
                TB_temp = TB_temp + BR(i,k,v_temp(i,k));
                MSE_over_temp = MSE_over_temp + MSE(i,k,v_temp(i,k))*N_ft(i,k)/sum(N_ft);
                if((F_vi(No_tile_vi) == i) &&(T_vi(No_tile_vi) == k))
                    Sum_BR_vi_temp = Sum_BR_vi_temp + BR(i,k,v_vi);
                else
                    Sum_BR_in_temp = Sum_BR_in_temp + BR(i,k,v_in);
                end
            end
        end
        if Sum_BR_in_temp == 0
            Sum_BR_in_temp = 1000000;
        end
        if((Sum_BR_vi_temp/No_tile_vi)/(Sum_BR_in_temp/No_tile_in) < rb) ...% Ratio of average bitrate
            && (TB_temp <= T_e) && (v_vi >=  v_sel_vi)
            Aver_BR_vi = Sum_BR_vi_temp/No_tile_vi;
            Aver_BR_in =Sum_BR_in_temp/No_tile_in;
            v_sel = v_temp;
            v_sel_vi = v_vi;
            TB_sel = TB_temp;
            MSE_over_sel = MSE_over_temp;
        end
    end
end

