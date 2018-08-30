%T(i) i=1:No_tile is selected versions for tile i
function [T,current_BW,m_,n_,P,B] = F_BellLab...
    (Fh, Fv, vp_W,vp_H, erp_W, erp_H, phi,theta,No_tile,No_ver,Bw,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H,Uti,Cti)

P=0* ones(1,No_tile);  % probality of tile include viewport = 0 or 1            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Viewport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m_,n_] = F_ExtractERPCodOfVP(Fh, Fv, vp_W,vp_H, erp_W, erp_H, phi,theta);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign value for P
for m=1:vp_W
    for n=1:vp_H
        for t=1:No_tile %tile
            if((n_(n,m)<HB_tile_H(t))&&(n_(n,m)>=LB_tile_H(t))...
                    &&(m_(n,m)<HB_tile_W(t))&&(m_(n,m)>=LB_tile_W(t)))
               P(t) = 1;
               break; % xem lai
            end
        end
    end
end

%% Find versions for each tile
% Sort (i,q) pairs in decreasing order of PU/C 
A = zeros(4, No_tile * No_ver);
for i=1:No_tile
    for q=1:No_ver
        A(1, (i-1)*No_ver + q) = i;
        A(2, (i-1)*No_ver + q) = q;
        A(3, (i-1)*No_ver + q) = (Uti(i,q) * P(i))/Cti(i,q);
        A(4, (i-1)*No_ver + q) = Cti(i,q);
    end
end
%
[ret, idx] = sort(A(3,:), 'descend');
tmp = A(1,:);
B(1,:) = tmp(idx);
tmp = A(2,:);
B(2,:) = tmp(idx);
tmp = A(3,:);
B(3,:) = tmp(idx);
tmp = A(4,:);
B(4,:) = tmp(idx);

%disp(A);
%disp(B);

current_BW = 0;
T = zeros(1,No_tile);
Min_Bw = 0;
for i=1:No_tile
	Min_Bw = Min_Bw + Cti(i,1);	
end
if Min_Bw <= Bw
	T = ones(1,No_tile);
	current_BW = Min_Bw;
end

%% Calculate the selected version for each tile
for j=1:No_tile * No_ver
    %fprintf('%d %d %d %d\n',B(1,j), B(2,j), B(4,j), current_BW);
    if current_BW + B(4,j) <= Bw && T(B(1,j)) == B(2,j) - 1
        T(B(1,j)) = B(2,j);
        current_BW = current_BW + B(4,j);
    end
end


 