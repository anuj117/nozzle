%% program to calculte grid points
function [X,Y,AR,max_t,yp]=MAIN(Me,n,gama,D)
nodes=n*(n+3)/2;
%preallocation of paramters
x=zeros(n,n+1);
y=zeros(n,n+1);
angle_mew=zeros(1,nodes);
angle_t_mat=zeros(n,n+1);
angle_mew_mat=zeros(n,n+1);
%pradtl meyer function and inverse prandtl meyer function
PM=@(M) (sqrt((gama+1)/(gama-1))*atan(sqrt((gama-1)*(M*M-1)/(gama+1))))-atan(sqrt(M*M-1));
function v=inv_PM(v) %nested function for inverse prandtl meyer function
        A=1.3604;
        B=0.0962;
        C=-0.5127;
        F=-0.6722;
        E=-0.3278;
        Vo=(pi/2)*((6).^(1/2)-1);
        y=(v/Vo).^(2/3);
        v=(1+A*y+B*y^2+C*y^3)/(1+F*y+E*y^2);       
 end
PM_inv=@(v) inv_PM(v);%function handle to inv_PM
max_wall=(PM(Me))*(180/(2*pi));%max wall angle 
d_angle=max_wall/n;%division of maxwall angle among each characteristic
[angle_t,v,KL,KR]=programx(max_wall,d_angle,n);

% KL is left running characteristics
% KR is right running charcteristics
%angle_t is angle theta for ach nodes
%% angle_mew prandtl meyer angle
for ii=1:nodes
    angle_mew(ii)=asind(1/PM_inv(v(ii)*pi/180));
end
ind=1;
for ii=1:n
    for kk=ii:n+1
        angle_t_mat(ii,kk)=angle_t(ind);
        angle_mew_mat(ii,kk)=angle_mew(ind);
        ind=ind+1;
    end
end
ind=1;
cp=0;
cl=0;
%%calculation of x,y points
for ii=1:n
    for kk=ii:n+1
        if ii==1%first right running charcteristic
            if kk==1%first point on first right
                x(ii,kk)=-D/(tand(angle_t_mat(ii,kk)-angle_mew_mat(ii,kk)));
                y=0;
            elseif (kk==n+1)%last wall point
                cp=(tand(0.5*(angle_t_mat(ii,kk-1)+angle_mew_mat(ii,kk-1)+angle_t_mat(ii,kk)+angle_mew_mat(ii,kk))));
                x(ii,kk)=(-D-x(ii,kk-1)*cp+y(ii,kk-1))/(tand(0.5*(max_wall+angle_t_mat(ii,kk)))-cp);
                y(ii,kk)=D+x(ii,kk)*tand(0.5*(max_wall+angle_t_mat(ii,kk)));
                
            else %mid points between wall and axis
                cp=(tand((0.5*(angle_t_mat(ii,kk-1)+angle_mew_mat(ii,kk-1)+angle_t_mat(ii,kk)+angle_mew_mat(ii,kk)))));
                x(ii,kk)=(-D+y(ii,kk-1)-cp*x(ii,kk-1))/(tand(angle_t_mat(ii,kk)-angle_mew_mat(ii,kk))-cp);
                y(ii,kk)=D+x(ii,kk)*tand(angle_t_mat(ii,kk)-angle_mew_mat(ii,kk));
                
            end
        
        %for rest of the characteristics
        else
            if(kk==ii)
                x(ii,kk)=x(ii-1,kk)-y(ii-1,kk)/(tand(0.5*(angle_t_mat(ii-1,kk)-angle_mew_mat(ii-1,kk)+angle_t_mat(ii,kk)-angle_mew_mat(ii,kk))));
                y(ii,kk)=0;
            elseif kk==n+1%points on wall
                mr=(tand(0.5*(angle_t_mat(ii,kk-1)+angle_mew_mat(ii,kk)+angle_t_mat(ii,kk-1)+angle_mew_mat(ii,kk))));
                ml=(tand(0.5*(angle_t_mat(ii-1,kk)+angle_t_mat(ii,kk))));
                x(ii,kk)=(y(ii-1,kk)-ml*x(ii-1,kk)-y(ii,kk-1)+mr*x(ii,kk-1))/(mr-ml);
                y(ii,kk)=y(ii-1,kk)+ml*(x(ii,kk)-x(ii-1,kk));
    
            else % mid points
                mr=(tand(0.5*(angle_t_mat(ii,kk-1)+angle_t_mat(ii,kk)+angle_mew_mat(ii,kk-1)+angle_mew_mat(ii,kk))));
                ml=(tand(0.5*(angle_t_mat(ii-1,kk)-angle_mew_mat(ii-1,kk)+angle_t_mat(ii,kk)-angle_mew_mat(ii,kk))));
                x(ii,kk)=(ml*x(ii-1,kk)+y(ii,kk-1)-mr*x(ii,kk-1)-y(ii-1,kk))/(ml-mr);
                y(ii,kk)=y(ii,kk-1)+mr*(x(ii,kk)-x(ii,kk-1));
                
            end
        end
    end
end
figure
hold on;
title('plot of all nodal points')
 for ii=1:n
       for kk=ii:n+1
          if ii==1
                plot([0,x(ii,kk)],[D,y(ii,kk)]);
          else
            plot([x(ii-1,kk),x(ii,kk)],[y(ii-1,kk),y(ii,kk)]);
          end
       end
 end
    for ii=1:n
        for kk=ii:n
            plot([x(ii,kk),x(ii,kk+1)],[y(ii,kk),y(ii,kk+1)]);
        end
    end
hold off;
X=zeros(1,n+1);
Y=zeros(1,n+1);
filename=fopen('points.dat','w');
ind=1;
X(1)=0;
Y(1)=D;
for ii=1:n
    X(ind+1)=x(ii,n+1);
    Y(ind+1)=y(ii,n+1);
    ind=ind+1;
end
figure 
hold on;
for i=1:n
    plot([X(i),X(i+1)],[Y(i),Y(i+1)]);
end
hold off;
for ii=1:n+1
    mat=[X(ii),Y(ii)];
    fprintf(filename,'%f %f\n',mat);
end
AR=Y(n+1)^2/D^2;
max_t=max_wall;
yp=Y(n+1);
fclose(filename);
end

