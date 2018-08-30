%% inputs:
% - w: width
% - h: height
% - view1: viewport 1 in yuv file
% - view1_ori: viewport 1 original in yuv file
% - view2: viewport 2 in yuv file
% - view2_ori: viewport 2 original in yuv file
% - x: the width in pixels of part from view 1
% - y: the width in pixels of part from view 2

%% output: 
% PSNR of the new viewport created from a part of x pixels of view
% 1 and a part of y pixels of view 2.

function [psnr, v1,v2,v1_ori,v2_ori]=F_viewportPSNR(w, h, view1, view1_ori, view2, view2_ori)
    %
    FRAME_NUM = 30; % number of video frames.
    START_FR_ID = 1;
    NUM_FR_READ = 30;
    M = 15;
    N = w/M;
    
    
    % load inputs 
    v1 = F_loadFileYuv(view1, w, h, FRAME_NUM, START_FR_ID, NUM_FR_READ, 0);
    v1_ori = F_loadFileYuv(view1_ori, w, h, FRAME_NUM, START_FR_ID, NUM_FR_READ, 0);
    v2 = F_loadFileYuv(view2, w, h, FRAME_NUM, START_FR_ID, NUM_FR_READ, 0);
    v2_ori = F_loadFileYuv(view2_ori, w, h, FRAME_NUM, START_FR_ID, NUM_FR_READ, 0);
    %     
    %     disp(v1);
    %     disp(v1_ori);
    %     disp(v2);
    %     disp(v2_ori);
    fname = sprintf('psnr_M=%d.txt',M);
    fout = fopen(fname, 'w');

    for i=1:NUM_FR_READ
        v1_tmp = v1(:,:,i);
        v1_ori_tmp = v1_ori(:,:,i);
        v2_tmp = v2(:,:,i);
        v2_ori_tmp = v2_ori(:,:,i);
        
        for j=1:N+1
            x = M * (j - 1);
            y = w - x;
            fprintf('%d:%d\n',x,y);
            % create view
            view = [v1_tmp(:,w-x+1:w) v2_tmp(:,1:y)];
            view_ori = [v1_ori_tmp(:,w-x+1:w) v2_ori_tmp(:,1:y)];

            % compute PSNR
            psnr(i,j) = 10*log10(255*255/immse(view, view_ori));
            fprintf(fout,'%.2f\t',psnr(i,j));
        end
        fprintf(fout, '\n');
    end
    