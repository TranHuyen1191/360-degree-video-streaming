%T(i) i=1:No_tile is selected versions for tile i
function [PSNR_error] = F_France_Cube_Error(Fh, Fv, vp_W,vp_H, ...
    face_W, face_H, phi,theta,No_face,tile_hori_num,tile_ver_num,...
    No_ver,v_sel,MSE,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MSE_over = 0;
    [N_ft,N_f,f_,m_,n_,m_c,n_c] = F_ExtractCubeTileCodOfVP(Fh, Fv,...
    vp_W,vp_H, face_W, face_H, phi,theta,tile_hori_num,...
    tile_ver_num,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H);
    for i = 1:No_face
        for k = 1:tile_hori_num*tile_ver_num
            MSE_over = MSE_over + MSE(i,k,v_sel(i,k))*N_ft(i,k)/sum(N_ft);
            PSNR_error = 10*log10(255*255/MSE_over);
        end
    end
end

