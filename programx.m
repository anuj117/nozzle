% function to calculate theta,kp,km,v for a given moc problems
function [t_angle,mew_angle,kl,kr]=programx(t_max,t_min,n)
delta_t=(t_max-t_min)/(n-1);
%% preallocation of varriables
nodes=n*(n+3)/2;
t_angle=zeros(1,nodes);
mew_angle=zeros(1,nodes);
kl=zeros(1,nodes);
kr=zeros(1,nodes);
%%calculating values at the nodes of first left running characteristics
for ii=1:n
    t_angle(ii)=t_min+(ii-1)*delta_t;
    mew_angle(ii)=t_angle(ii);
    kr(ii)=t_angle(ii)+mew_angle(ii);
    kl(ii)=t_angle(ii)-mew_angle(ii);
end
cur_ind=n+1;
% % calculating values at first point of the contour wall from throat
t_angle(cur_ind)=t_angle(n);
mew_angle(cur_ind)=t_angle(n);
kr(cur_ind)=kr(n);
kl(cur_ind)=kl(n);
cur_ind=cur_ind+1;
%%for remaining characteristics after the allocation of values at nodes of
%%first right running characteristic
for jj=1:n-1 
    t_angle(cur_ind)=0;
    kr(cur_ind)=kr(1+jj);
    mew_angle(cur_ind)=kr(cur_ind)-t_angle(cur_ind);
    kl_all=t_angle(cur_ind)-mew_angle(cur_ind);
    kl(cur_ind)=kl_all;
    cur_ind=cur_ind+1;
    for nodes=1:n-jj-1
        kr(cur_ind)=kr(jj+1+nodes);
        kl(cur_ind)=kl_all;
        t_angle(cur_ind)=(kr(cur_ind)+kl(cur_ind))/2;
        mew_angle(cur_ind)=(kr(cur_ind)-kl(cur_ind))/2;
        zx=kr(cur_ind)+kl(cur_ind);
        cur_ind=cur_ind+1;
        
    end
    %for nodes on nozzle wall
    kr(cur_ind)=kr(cur_ind-1);
    t_angle(cur_ind)=t_angle(cur_ind-1);
    mew_angle(cur_ind)=mew_angle(cur_ind-1);
    kr(cur_ind)=kr(cur_ind-1);
    kl(cur_ind)=kl(cur_ind-1);
    cur_ind=cur_ind+1;
end
end