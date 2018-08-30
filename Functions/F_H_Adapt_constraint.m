%metric = xichma_of_all_tile(delta MSE * No_visiblepixel / deltabirate) / No_pixel_of_vp
%v_sel(i) i=1:No_face    are selected versions for face i
function [TD_sel,v_sel,TB_sel,N_ft] = F_H_Adapt_constraint(Fh, Fv, vp_W,vp_H, face_W, ...
    face_H, phi,theta,No_face,tile_hori_num,tile_ver_num,No_ver,T_e,BR,...
    MSE,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Viewport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N_ft] = F_ExtractCubeTileCodOfVP...
    (Fh, Fv, vp_W,vp_H, face_W, face_H, phi,theta,tile_hori_num,...
    tile_ver_num,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% Calculate the selected version for each tile
v_sel = ones(No_face,tile_hori_num*tile_ver_num); % array of selected versions
v_temp = ones(No_face,tile_hori_num*tile_ver_num); % initialize 1 (no transmit)
TB_sel =0;
TB_temp =0;
TD_sel =65025;
TD_temp =0;
TD_sel_min = 0;
Min_T_e = 0; % T_e min to transmit all face with lowest versions
N_ft_max = 0; % the highest weight
f_Nmax = 0; % Face having the highest weight
t_Nmax = 0; % Tile having the highest weight

%% if T_e is enough --> transmit all face
for i = 1:No_face
	for k = 1:tile_hori_num*tile_ver_num
        % find the highest weight to detect G1 and G2
        if(N_ft(i,k) > N_ft_max)
            N_ft_max = N_ft(i,k);
            f_Nmax = i;
            t_Nmax = k;
        end
		Min_T_e = Min_T_e + BR(i,k,2);
        TD_sel_min = TD_sel_min + MSE(i,k,2) * N_ft(i,k)/(vp_W*vp_H);
	end
end

if T_e >= Min_T_e 
	v_sel = ones(No_face,tile_hori_num*tile_ver_num)*2;
	TB_sel = Min_T_e;
    TD_sel = TD_sel_min;
end

%% Detect G1 and G2
No_g1 = 0; % Number of tiles in G1
No_g2 = 0; % Number of tiles in G2
for i = 1:No_face
	for k = 1:tile_hori_num*tile_ver_num 
        if(N_ft(i,k) > N_ft_max/2)
            No_g1 = No_g1 +1;
            Gf_1(No_g1) = i;
            Gt_1(No_g1) = k;
        elseif (N_ft(i,k) > 0)
            No_g2 = No_g2 +1;
            Gf_2(No_g2) = i;
            Gt_2(No_g2) = k;
        end
	end
end

%% find selected versions for all tiles
v_temp = v_sel;
for v2 = 2:No_ver+1
    if No_g2  > 0 
        for cnt_g2 = 1:No_g2
            v_temp(Gf_2(cnt_g2),Gt_2(cnt_g2)) = v2;
        end
    end
    for v1 = v2 : No_ver+1
        for cnt_g1 = 1:No_g1
            v_temp(Gf_1(cnt_g1),Gt_1(cnt_g1)) = v1;
        end
        TB_temp  = 0;
        TD_temp = 0;
        % Calculate TB and TD
        for i = 1:No_face
            for k = 1:tile_hori_num*tile_ver_num 
                TB_temp = TB_temp + BR(i,k,v_temp(i,k));
                TD_temp = TD_temp + MSE(i,k,v_temp(i,k)) * N_ft(i,k)/(vp_W*vp_H);
            end
        end
        if(TD_temp <= TD_sel) && (TB_temp <= T_e)
			% assign versions of tiles 
            v_sel = v_temp;
            TB_sel = TB_temp;
            TD_sel = TD_temp;
        end
    end
end


%% Process the remaining bandwidth
v_Nmax = v_sel(f_Nmax,t_Nmax);
TB_woNmax = TB_sel - BR(f_Nmax,t_Nmax,v_Nmax);
TD_woNmax = TD_sel - MSE(f_Nmax,t_Nmax,v_Nmax)* N_ft(f_Nmax,t_Nmax)/(vp_W*vp_H);
for v_ = v_Nmax +1 : No_ver+1
    % Calculate TB and TD
    TB_temp = TB_woNmax + BR(f_Nmax,t_Nmax,v_);
    TD_temp = TD_woNmax + MSE(f_Nmax,t_Nmax,v_)* N_ft(f_Nmax,t_Nmax)/(vp_W*vp_H);
    if( TB_temp <= T_e)
        TB_sel = TB_temp;
        TD_sel = TD_temp;
        v_sel(f_Nmax,t_Nmax) = v_ ;
    end
end

for i = 1:No_face
    for k = 1:tile_hori_num*tile_ver_num
        if(N_ft(i,k) > 0)
            if v_sel(i,k)+1 <= No_ver +1
                TB_temp = TB_sel - BR(i,k,v_sel(i,k)) + BR(i,k,v_sel(i,k) +1);
                TD_temp = TD_sel - (MSE(i,k,v_sel(i,k)) + MSE(i,k,v_sel(i,k) +1))* N_ft(i,k)/(vp_W*vp_H);
                if (TB_temp <= T_e)
                    v_sel(i,k) = v_sel(i,k)+1;
                    TB_sel = TB_temp;
                    TD_sel = TD_temp;
                end
            end
        end
    end
end

toc
