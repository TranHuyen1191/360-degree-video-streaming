function [N_ft,N_f,f_,m_f,n_f,m_c,n_c] = F_ExtractCubeTileCodOfVP(Fh, Fv, vp_W,vp_H, face_W, face_H, phi,theta,tile_hori_num,tile_ver_num,LB_tile_W,LB_tile_H,HB_tile_W,HB_tile_H)
% f = 1:6, m_f = 0:face_W-1, m_c = 0: face_W*4-1
% Using to extract cordinates in Cubemap coressponding to points on Viewport.
%Point (i,j) on viewport coresponding to (m_f(i,j),n_f(i,j)) on face f_ in Cubemap image.
%Point (i,j) on viewport coresponding to (m_c(i,j),n_c(i,j)) on Cubemap image 4x3.
% i and m_f are width 
% j and n_f are height
% N_f(f): number of pixels on viewport rounded from face f.
% N_ft(f,t): number of pixels on viewport rounded from tile t of face f.
m_c=zeros(vp_H,vp_W);
n_c=zeros(vp_H,vp_W);
m_f=zeros(vp_H,vp_W);
n_f=zeros(vp_H,vp_W);
f_=zeros(vp_H,vp_W);
N_f =zeros(1,6);
N_ft =zeros(6,tile_hori_num*tile_ver_num);
%   fname = sprintf('out_%.0f_%.0f.txt' ,phi* 180/pi,theta* 180/pi);
%   fout = fopen(fname,'w');
%   fprintf(fout, 'phi  theta	m	n	phi_    theta_	X	Y	Z	x_	y_	z   m   n   \n');
%fname_log = sprintf('Log_positionVP_facecube_%f.txt',phi);
%fout_log = fopen(fname_log,'w');
%fprintf(fout_log,'position on VP --> postion on cube\n');
%fprintf(fout_log,'m n face m n\n');

for m=0:(vp_W-1) % width position of point on viewports
    for n=0:(vp_H-1) % height position of point on viewports

        u = (m+0.5) * 2 * tan(Fh/2)/vp_W;
        v = (n+0.5) * 2 * tan(Fv/2)/vp_H;

        x = u - tan(Fh/2);
        y = -v + tan(Fv/2);
        z = 1;

        x_ = x/sqrt(x*x + y*y + z*z);
        y_ = y/sqrt(x*x + y*y + z*z);
        z_ = z/sqrt(x*x + y*y + z*z);

        R = [cos(phi + pi/2)   -sin(phi + pi/2)*sin(theta)     sin(phi+pi/2) * cos(theta);
             0                 cos(theta)                      sin(theta)                ;
             -sin(phi + pi/2)  -cos(phi + pi/2) * sin(theta)   cos(phi + pi/2) * cos(theta)];

        A= R*[x_;y_;z_];
        X=A(1);
        Y=A(2);
        Z=A(3);

        % cubemap
        if (abs(X)>= abs(Y)) &&  (abs(X)>= abs(Z)) &&  ((X)>0) 
            f = 0;
            u_ = -Z/abs(X);
            v_ = -Y/abs(X);
        elseif (abs(X)>= abs(Y)) &&  (abs(X)>= abs(Z)) &&  ((X)<0) 	
            f = 1;
            u_ = Z/abs(X);
            v_ = -Y/abs(X);
        elseif (abs(Y)>= abs(X)) &&  (abs(Y)>= abs(Z)) &&  ((Y)>0) 	
            f = 2;
            u_ = X/abs(Y);
            v_ = Z/abs(Y);
        elseif (abs(Y)>= abs(X)) &&  (abs(Y)>= abs(Z)) &&  ((Y)<0) 	
            f = 3;
            u_ = X/abs(Y);
            v_ = -Z/abs(Y);
        elseif (abs(Z)>= abs(X)) &&  (abs(Z)>= abs(Y)) &&  ((Z)>0) 	
            f = 4;
            u_ = X/abs(Z);
            v_ = -Y/abs(Z);
        elseif (abs(Z)>= abs(X)) &&  (abs(Z)>= abs(Y)) &&  ((Z)<0) 	
            f = 5;
            u_ = -X/abs(Z);
            v_ = -Y/abs(Z);
        end
        m_f(n+1,m+1) = int32((u_+1) * face_W/2 - 0.5);
        n_f(n+1,m+1) = int32((v_+1) * face_H/2 - 0.5);
        f_(n+1,m+1) = f+1;
        N_f(f+1) = N_f(f+1)+1;
        %if ((f_(n+1,m+1) == 3))  
        %    fprintf(fout_log,'%d %d %d %f %f \n',m,n,f_(n+1,m+1), m_f(n+1,m+1),n_f(n+1,m+1));
        %end
        % calculate number of pixels for each tile
		for t=1:tile_hori_num*tile_ver_num %tile
			if((n_f(n+1,m+1)<HB_tile_H(t))&&(n_f(n+1,m+1)>=LB_tile_H(t))&&(m_f(n+1,m+1)<HB_tile_W(t))&&(m_f(n+1,m+1)>=LB_tile_W(t)))
			   N_ft(f+1,t) = N_ft(f+1,t) + 1;
			   break; % xem lai
			end
        end
        if f == 0
           m_c(n+1,m+1) = m_f(n+1,m+1)+face_W*2;
           n_c(n+1,m+1) = n_f(n+1,m+1)+face_H; 
        elseif f == 1
           m_c(n+1,m+1) = m_f(n+1,m+1);
           n_c(n+1,m+1) = n_f(n+1,m+1)+face_H;
        elseif f == 2
           m_c(n+1,m+1) = m_f(n+1,m+1)+face_W;
           n_c(n+1,m+1) = n_f(n+1,m+1);
        elseif f == 3
           m_c(n+1,m+1) = m_f(n+1,m+1)+face_W;
           n_c(n+1,m+1) = n_f(n+1,m+1)+face_H*2;
        elseif f == 4
           m_c(n+1,m+1) = m_f(n+1,m+1)+face_W;
           n_c(n+1,m+1) = n_f(n+1,m+1)+face_H;
        elseif f == 5
           m_c(n+1,m+1) = m_f(n+1,m+1)+face_W*3;
           n_c(n+1,m+1) = n_f(n+1,m+1)+face_H;
        end    
    end
end

%   fclose(fout);

