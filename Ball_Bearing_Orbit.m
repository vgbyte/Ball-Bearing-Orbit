clear all; 

% INPUT DATA ( All the dimensions are in SI unit )

    Wi  = 180; Wo = 0;      % Inner and Outer race angular velocity in rad/sec, Wi=5rad/sec
    W   = Wi/2;           % Angular velocity of ball cenetr rotating about bearing axis 
    Rb  = 3.967e-3;       % Ball Radius in m
    Ri  = 13.281e-3;      % Inner race radius in m
    Ro  = 21.226e-3;      % Outer race radius in m
    N   = 8;              % Number of balls
    Mb  = 0.002;          % Single Ball Mass
    Mi  = 0.045;          % Inner race mass
    Mo  = 0.05;           % Outer race mass
    c_s = 740;            % Damping coefficient in N-s/m
    c_b = 0;%800;            % Damping coefficient in N-s/m
    Kpb = 2.8397e+05;     % Hertzian Constant Coefficient in N/m^(3/2)
    K_s = 3.7e7;          % Linear stiffness between shaft and bearing in N/m
    
    ts  = 1e-6;           % Time step in sec
    TT  = 3;              % Total time in sec
    Fo  = 10;             % Harmonic force magnitude in N
    
    maxItr=TT/ts;         % Total iteration 
    
    x_b=zeros(N,maxItr+1);
    y_b=zeros(N,maxItr+1);
    dr_b=zeros(N,maxItr+1);
    x_i=zeros(1,maxItr+1);
    x_o=zeros(1,maxItr+1);
    y_i=zeros(1,maxItr+1);
    y_o=zeros(1,maxItr+1);
    u_i=zeros(1,maxItr+1);
    u_o=zeros(1,maxItr+1);
    v_i=zeros(1,maxItr+1);
    v_o=zeros(1,maxItr+1);
    t  =zeros(1,maxItr+1);
    
    X=zeros(1,maxItr);
    Y=zeros(1,maxItr);
    
    
% Initial position & velocity vector of ball and races center
    x_i(1)=0; y_i(1)=0;
    u_i(1)=0; v_i(1)=0;
    x_o(1)=0; y_o(1)=0;
    u_o(1)=0; v_o(1)=0;
    for n=1:N
        % Position vector of ball center
        x_b(n,1) = (Ri+Rb)*cos((2*pi/N)*(n-1));
        y_b(n,1) = (Ri+Rb)*sin((2*pi/N)*(n-1));
        dr_b(n,1)= 0;      % Only radial velocity is defined
    end
    

for i=1:maxItr
    F_inner=[0,0];
    F_outer=[0,0];
    for n=1:N
        F_damping_ball_inner = 0; F_spring_ball_inner = 0;
        F_damping_ball_outer = 0; F_spring_ball_outer = 0;
        F_spring_ball_mag    = 0; F_damping_ball_mag  = 0;
        r=((x_i(i)-x_b(n,i))^2+(y_i(i)-y_b(n,i))^2)^(0.5);
        dr=dr_b(n,i);
        % Check whether their is inner race and ball deformation
        del_i=Ri+Rb-r;
        if del_i>0
            F_spring_ball_inner=Kpb*del_i^(1.5);
            r_b=[x_b(n,i),y_b(n,i)];
            r_c1=[x_i(i),y_i(i)];
            v_c1=[u_i(i),v_i(i)];
            v_b=dr*(r_b-r_c1)/norm(r_b-r_c1);
            r_c1_b=r_b-r_c1;
            v_c1_b=v_b-v_c1;
            x_dot=projection_A_on_B(v_c1_b,r_c1_b);
            F_damping_ball_inner=-x_dot*c_b;
        end
        
        % Check whether their is outer race and ball deformation
        del_o=Rb-Ro+((x_o(i)-x_b(n,i))^2+(y_o(i)-y_b(n,i))^2)^(0.5);
        if del_o>0
            r_b=[x_b(n,i),y_b(n,i)];
            r_c1=[x_i(i),y_i(i)];
            r_c2=[x_o(i),y_o(i)];
            v_c2=[u_o(i),v_o(i)];
            v_b=dr*(r_b-r_c1)/norm(r_b-r_c1);
            r_c1_b=r_c1-r_b;
            r_c2_b=r_c2-r_b;
            v_c2_b=v_c2-v_b;
            
            x_dot=projection_A_on_B(v_c2_b,r_c2_b);
            F_damping_ball_mag=-x_dot*c_b;
            F_spring_ball_mag=-Kpb*del_o^(1.5);
            cos_theta=dot(r_c1_b,r_c2_b)/(norm(r_c1_b)*norm(r_c2_b));
            F_damping_ball_outer=F_damping_ball_mag*cos_theta;
            F_spring_ball_outer=F_spring_ball_mag*cos_theta;
           
        end
        
        F_damping_ball_total=F_damping_ball_inner+F_damping_ball_outer;
        F_spring_ball_total=F_spring_ball_inner+F_spring_ball_outer;
        F_ball_total=F_damping_ball_total+F_spring_ball_total;
        
        % RK4 for each ball
        kb1_1 = ts*fun_distance(dr);
        kb2_1 = ts*fun_ball(Mb,W,F_ball_total,r);
        kb1_2 = ts*fun_distance(dr+0.5*kb2_1);
        kb2_2 = ts*fun_ball(Mb,W,F_ball_total,r+0.5*kb1_1);
        kb1_3 = ts*fun_distance(dr+0.5*kb2_2);
        kb2_3 = ts*fun_ball(Mb,W,F_ball_total,r+0.5*kb1_2);
        kb1_4 = ts*fun_distance(dr+kb2_3);
        kb2_4 = ts*fun_ball(Mb,W,F_ball_total,r+kb1_3);
        
        r_new = r +(kb1_1 + 2*kb1_2 + 2*kb1_3 + kb1_4)/6;
        dr_new= dr+(kb2_1 + 2*kb2_2 + 2*kb2_3 + kb2_4)/6;
        
        theta=atan2(y_b(n,i)-y_i(i),x_b(n,i)-x_i(i))+2*pi;
        x_b(n,i+1)=r_new*cos(theta+W*ts);
        y_b(n,i+1)=r_new*sin(theta+W*ts);
        dr_b(n,i+1)=dr_new;
        
        if n==1
            X(i)=x_b(n,i);
            Y(i)=y_b(n,i);
        end
        
        r_b=[x_b(n,i),y_b(n,i)];
        r_c1=[x_i(i),y_i(i)];
        r_c2=[x_o(i),y_o(i)];
        r_c1_b=r_c1-r_b;
        r_c2_b=r_c2-r_b;
        F_inner=F_inner+(F_spring_ball_inner+F_damping_ball_inner)*r_c1_b/norm(r_c1_b); 
        F_outer=F_outer+(F_spring_ball_mag+F_damping_ball_mag)*r_c2_b/norm(r_c2_b);
    end
    
    
    kix1_1 = ts*fun_distance(u_i(i));
    kix2_1 = ts*fun_race_x(Mi,c_s,K_s,F_inner(1),Fo,Wi,u_i(i),x_i(i),t(i));
    kix1_2 = ts*fun_distance(u_i(i)+0.5*kix2_1);
    kix2_2 = ts*fun_race_x(Mi,c_s,K_s,F_inner(1),Fo,Wi,u_i(i)+0.5*kix2_1,x_i(i)+0.5*kix1_1,t(i)+ts/2);
    kix1_3 = ts*fun_distance(u_i(i)+0.5*kix2_2);
    kix2_3 = ts*fun_race_x(Mi,c_s,K_s,F_inner(1),Fo,Wi,u_i(i)+0.5*kix2_2,x_i(i)+0.5*kix1_2,t(i)+ts/2);
    kix1_4 = ts*fun_distance(u_i(i)+kix2_3);
    kix2_4 = ts*fun_race_x(Mi,c_s,K_s,F_inner(1),Fo,Wi,u_i(i)+kix2_3,x_i(i)+kix1_3,t(i)+ts);
    
    kiy1_1 = ts*fun_distance(v_i(i));
    kiy2_1 = ts*fun_race_y(Mi,c_s,K_s,F_inner(2),Fo,Wi,v_i(i),y_i(i),t(i));
    kiy1_2 = ts*fun_distance(v_i(i)+0.5*kiy2_1);
    kiy2_2 = ts*fun_race_y(Mi,c_s,K_s,F_inner(2),Fo,Wi,v_i(i)+0.5*kiy2_1,y_i(i)+0.5*kiy1_1,t(i)+ts/2);
    kiy1_3 = ts*fun_distance(v_i(i)+0.5*kiy2_2);
    kiy2_3 = ts*fun_race_y(Mi,c_s,K_s,F_inner(2),Fo,Wi,v_i(i)+0.5*kiy2_2,y_i(i)+0.5*kiy1_2,t(i)+ts/2);
    kiy1_4 = ts*fun_distance(v_i(i)+kiy2_3);
    kiy2_4 = ts*fun_race_y(Mi,c_s,K_s,F_inner(2),Fo,Wi,v_i(i)+kiy2_3,y_i(i)+kiy1_3,t(i)+ts);
    
    kox1_1 = ts*fun_distance(u_o(i));
    kox2_1 = ts*fun_race_outer(Mo,F_outer(1));
    kox1_2 = ts*fun_distance(u_o(i)+0.5*kox2_1);
    kox2_2 = ts*fun_race_outer(Mo,F_outer(1));
    kox1_3 = ts*fun_distance(u_o(i)+0.5*kox2_2);
    kox2_3 = ts*fun_race_outer(Mo,F_outer(1));
    kox1_4 = ts*fun_distance(u_o(i)+kox2_3);
    kox2_4 = ts*fun_race_outer(Mo,F_outer(1));
    
    koy1_1 = ts*fun_distance(v_o(i));
    koy2_1 = ts*fun_race_outer(Mo,F_outer(2));
    koy1_2 = ts*fun_distance(v_o(i)+0.5*koy2_1);
    koy2_2 = ts*fun_race_outer(Mo,F_outer(2));
    koy1_3 = ts*fun_distance(v_o(i)+0.5*koy2_2);
    koy2_3 = ts*fun_race_outer(Mo,F_outer(2));
    koy1_4 = ts*fun_distance(v_o(i)+koy2_3);
    koy2_4 = ts*fun_race_outer(Mo,F_outer(2));
    
    
    % Updating position and velocity vector of inner race center
    x_i(i+1) = x_i(i)+(kix1_1 + 2*kix1_2 + 2*kix1_3 + kix1_4)/6;
    u_i(i+1) = u_i(i)+(kix2_1 + 2*kix2_2 + 2*kix2_3 + kix2_4)/6;
    y_i(i+1) = y_i(i)+(kiy1_1 + 2*kiy1_2 + 2*kiy1_3 + kiy1_4)/6;    
    v_i(i+1) = v_i(i)+(kiy2_1 + 2*kiy2_2 + 2*kiy2_3 + kiy2_4)/6;
    
    % Updating position and velocity vector of inner race center
    x_o(i+1) = x_o(i)+(kox1_1 + 2*kox1_2 + 2*kox1_3 + kox1_4)/6;
    u_o(i+1) = u_o(i)+(kox2_1 + 2*kox2_2 + 2*kox2_3 + kox2_4)/6;
    y_o(i+1) = y_o(i)+(koy1_1 + 2*koy1_2 + 2*koy1_3 + koy1_4)/6;    
    v_o(i+1) = v_o(i)+(koy2_1 + 2*koy2_2 + 2*koy2_3 + koy2_4)/6;
    
    t(i+1)=t(i)+ts;
end
figure(1)
plot(x_o,y_o)
title('Orbit of center of outer race')
grid on

figure(2)
plot(x_i,y_i)
title('Orbit of center of inner race')
grid on

figure(3)
plot(X,Y)
title('Path of center of Ball')
grid on

function a = fun_distance(v)
    a = v;
end

function a = fun_ball(Mb,W,F_ball_total,r)
    a=F_ball_total/Mb+r*W^2;
end

function a= fun_race_outer(M,F)
    a=F/M;
end

function a = fun_race_x(M,c,k,F_total,Fo,Wi,u,x,t)
    a=(F_total+Fo*cos(Wi*t)-c*u-k*x)/M;
end

function a = fun_race_y(M,c,k,F_total,Fo,Wi,v,y,t)
    a=(F_total+Fo*sin(Wi*t)-c*v-k*y)/M;
end

function a = projection_A_on_B(A,B)
    a=dot(A,B)/norm(B);
end