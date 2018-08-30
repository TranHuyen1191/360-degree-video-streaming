function [m_,n_] = F_ExtractERPCodOfVP(Fh, Fv, vp_W,vp_H, erp_W, erp_H, phi,theta)
% Using to extract cordinates in ERP coressponding to points on Viewport.
%Point (i,j) on viewport coresponding to (m_(i,j),n(i,j)) on ERP image.s
% i and m_ are width 
% j and n_ are height
m_=zeros(vp_H,vp_W);
n_=zeros(vp_H,vp_W);

%   fname = sprintf('out_%.0f_%.0f.txt' ,phi* 180/pi,theta* 180/pi);
%   fout = fopen(fname,'w');
%   fprintf(fout, 'phi  theta	m	n	phi_    theta_	X	Y	Z	x_	y_	z   m   n   \n');

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

        % erp
        phi_ = atan(-Z/X);
        theta_ = asin(Y/sqrt(X*X + Y*Y + Z*Z));
        if X < 0 % phi_  >180 degree
            phi_=pi+phi_;
        end

        
        u_ = phi_/(2*pi) + 0.5;
        v_ = 0.5 - theta_ / pi;
        
        m_temp =  int32(u_ * erp_W - 0.5);
        if m_temp > erp_W 
            m_(n+1,m+1) = m_temp-erp_W;
        else
            m_(n+1,m+1) = m_temp;
            
        end    
        n_(n+1,m+1) = int32(v_ * erp_H - 0.5);
        %   fprintf(fout, '%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n' ,phi* 180/pi,theta* 180/pi,m,n,phi_* 180/pi,theta_* 180/pi,X,Y,Z,x_,y_,z_, m_(n+1,m+1), n_(n+1,m+1));
    end
end
%   fclose(fout);

