imax=25;
jmax=25;
[u0,v0,p0]=Predictor_Corrector(imax,jmax,0.00002,0.1);

[u1,v1,p1]=Predictor_Corrector(imax,jmax,0.0002,1);

Velocity_Plots(imax,jmax,u0,v0,u1,v1)


function [u,v,p] = Predictor_Corrector(imax,jmax,dt,Re)
    % Mesh
    dx=1/(imax-1);
    dy=1/(jmax-1);
    x=0:dx:1.0;
    y=0:dy:1.0;

    % Time Simulation
    t=0.0;
    Fractional_tStep=0.02/dt;
    fprintf('Fr_dt=%d',Fractional_tStep);
    counter=1;

    Lid_Velocity=1;
    
    % Variable Declaration
    p=zeros(imax+1,jmax+1);
    u=zeros(imax,jmax+1);
    v=zeros(imax+1,jmax);

    % Laplacian Operator   
    [Laplacian_to_solve]=Laplacian(imax,jmax,dx,dy);    
    iter=0;
    while(t<=0.1)
    
        [us,vs]=velocities_predictor(imax,jmax,dt,dx,dy,Re,u,v,Lid_Velocity);    
        [R] = R_vector(imax,jmax,us,vs,dx,dy,dt);    
        
        [p]=Poisson(imax,jmax,Laplacian_to_solve,R);    
        
        [u,v]=velocities_corrector(imax,jmax,dt,dx,dy,us,vs,p,Lid_Velocity);    
        if iter-1==floor(Fractional_tStep*counter)
            counter=counter+1;
            Result_plots(x,y,imax,jmax,dx,v,p,t,Re)
        end
        iter=iter+1;
        t=t+dt;
    end
end


function [Laplacian_Matrix] = Laplacian(imax,jmax,dx,dy)
    nx=imax-1;
    ny=jmax-1;
    dx_i=1/dx;
    dy_i=1/dy;
    vec=[(2*dx_i^2)+(2*dy_i^2),-dx_i^2,zeros(1,nx-2),-dy_i^2,zeros(1,(nx*ny)-(3+nx-2))];
    Laplacian_Matrix=toeplitz(vec); 
    
    for j=1:ny
    	for i=1:nx
     
    		rw=i+(j-1)*nx;
    	
    		if(mod(rw,nx)==0) 
    	
    			x=nx;
    		else
    			x=mod(rw,nx);
    		end
     
    		if(mod(rw,ny)==0) 
    
    			y=floor(rw/ny);
    		else
    			y=floor(rw/ny)+1;
    		end
     
    		if(x==1||x==nx ) 
    			Laplacian_Matrix(rw,rw)=Laplacian_Matrix(rw,rw)-dx_i^2;
    		end
     
    		if(y==1 ||y==ny) 
    			Laplacian_Matrix(rw,rw)=Laplacian_Matrix(rw,rw)-dy_i^2;
    		end
     
    		if(mod(rw,nx)==1 &&rw~=1)
    			Laplacian_Matrix(rw,rw-1)=0;
    		end
    		if(mod(rw,nx)==0 && rw~=nx*ny) 
    			Laplacian_Matrix(rw,rw+1)=0;
    		end
    	end
    end
    Laplacian_Matrix(1,:)=0;
    Laplacian_Matrix(1,1)=1;
    
    return
end

% Predictor
function [ us,vs ] = velocities_predictor(imax,jmax,dt,dx,dy,Re,u,v,Lid_Velocity)
    us=zeros(imax,jmax+1);
    for i=2:imax-1
        for j=2:jmax
            us(i,j)=u(i,j)+dt*(-(((0.5*(u(i+1,j)+u(i,j)))^2-(0.5*(u(i-1,j)+u(i,j)))^2)/dx+(((0.5*(u(i,j+1)+u(i,j)))*(0.5*(v(i+1,j)+v(i,j))))-((0.5*(u(i,j-1)+u(i,j)))*(0.5*(v(i,j-1)+v(i+1,j-1)))))/dy)+((((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2))/Re));
        end
    end

    us(1,2:jmax)=0.0;
    us(imax,2:jmax)=0.0;
    us(1:imax,1)=-us(1:imax,2);
    us(1:imax,jmax+1)=2.0*Lid_Velocity-us(1:imax,jmax);

   
    vs=zeros(imax+1,jmax);
    for i=2:imax
        for j=2:jmax-1
            vs(i,j)=v(i,j)+dt*(-(((((0.5*(u(i,j+1)+u(i,j)))*(0.5*(v(i+1,j)+v(i,j))))-((0.5*(u(i-1,j+1)+u(i-1,j)))*(0.5*(v(i-1,j)+v(i,j)))))/dx)+((((0.5*(v(i,j+1)+v(i,j)))^2)-((0.5*(v(i,j-1)+v(i,j)))^2))/dy))+((((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)+((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2))/Re));
        end
    end

    vs(1,1:jmax)=-vs(2,1:jmax);
    vs(imax+1,1:jmax)=-vs(imax,1:jmax);
    vs(2:imax,1)=0.0;
    vs(2:imax,jmax)=0.0;
    return
end


function [ R ] = R_vector(imax,jmax,us,vs,dx,dy,dt)
    R_vec=zeros(imax+1,jmax+1);
    
    for i=2:imax
        for j=2:jmax
            R_vec(i,j)=-((us(i,j)-us(i-1,j))/dt/dx+(vs(i,j)-vs(i,j-1))/dt/dy);
        end
    end
   
    R=zeros((imax-1)*(jmax-1));
    index=1;
    
    for j=2:jmax
        for i=2:imax
            R(index)=R_vec(i,j);
            index=index+1;
        end
    end
    return
end


function [ p ] = Poisson(imax,jmax,Laplacian_to_solve,R)
    p=zeros(imax+1,jmax+1);
   
    pressure=Laplacian_to_solve\R;
    index=1;
    for j=2:jmax
        for i=2:imax
            p(i,j)=pressure(index);
            index=index+1;
        end
    end
    
   
    p(2:imax,1)=p(2:imax,2);
    p(2:imax,jmax+1)=p(2:imax,jmax);
    p(1,2:jmax)=p(2,2:jmax);
    p(imax+1,2:jmax)=p(imax,2:jmax);
    return
end


function [ u,v ] = velocities_corrector(imax,jmax,dt,dx,dy,us,vs,p,Lid_Velocity)
    u=zeros(imax,jmax+1);
    v=zeros(imax+1,jmax);
    
    % Corrector for u
    for i=2:imax-1
        for j=2:jmax
            u(i,j)=us(i,j)-(p(i+1,j)-p(i,j))*dt/dx;
        end
    end
    
    u(1,2:jmax)=0.0;
    u(imax,2:jmax)=0.0;
    u(2:imax-1,1)=-u(2:imax-1,2);
    u(2:imax-1,jmax+1)=2*Lid_Velocity-u(2:imax-1,jmax);
    
    % Corrector for v
    for i=2:imax
        for j=2:jmax-1
           v(i,j)=vs(i,j)-(p(i,j+1)-p(i,j))*dt/dy;
        end
    end
    
    v(2:imax,1)=0.0;
    v(2:imax,jmax)=0.0;
    v(1,2:jmax-1)=-v(2,2:jmax-1);
    v(imax+1,2:jmax-1)=-v(imax,2:jmax-1);
    return
end


function [] = Velocity_Plots(imax,jmax,u0,v0,u1,v1)
    dx=1/(imax-1);
    dy=1/(jmax-1);
    x=0:dx:1.0;
    y=0:dy:1.0;
    uc0=zeros(imax,jmax+1);
    vc0=zeros(imax+1,jmax);
    uc1=zeros(imax,jmax+1);
    vc1=zeros(imax+1,jmax);
    
    for j=1:jmax
        for i=1:imax
            uc0(i,j)=0.5*(u0(i,j+1)+u0(i,j));
            uc1(i,j)=0.5*(u1(i,j+1)+u1(i,j));
            vc0(i,j)=0.5*(v0(i+1,j)+v0(i,j));
            vc1(i,j)=0.5*(v1(i+1,j)+v1(i,j));
        end
    end
     
    x_center=zeros(length(x)-1,1);
    y_center=zeros(length(y)-1,1);
    for i=1:length(x_center)
        x_center(i)=(x(i)+x(i+1))/2;
    end
    for j=1:length(y_center)
        y_center(j)=(y(j)+y(j+1))/2;
    end
    u_center0=uc0(floor(imax/2),:);
    u_center1=uc1(floor(imax/2),:);
    
    figure
    plot(y_center,u_center0(1,2:end-1),'-go','MarkerSize',9,'DisplayName','Re=0.1');
    hold on;
    title('Normalized u comparision for Re=0.1 nd Re=1'); 
    xlabel('Y/L');
    ylabel('u/U');
    plot(y_center,u_center1(1,2:end-1),'-rx','DisplayName','Re=1');
    hold off;
    legend('Re=0.1','Re=1');
    saveas(gcf,'uc_velocity.jpg');
    v_center0=vc0(:,floor(jmax/2));
    v_center1=vc1(:,floor(jmax/2));
    
    figure
    plot(x_center',v_center0(2:end-1,1),'-go','MarkerSize',9,'DisplayName','Re=0.1');
    hold on;
    title('Normalized v comparision for Re=0.1 nd Re=1'); 
    xlabel('x/L');
    ylabel('v/U');
    plot(x_center',v_center1(2:end-1,1),'-rx','DisplayName','Re=1');
    hold off;
    legend('Re=0.1','Re=1');
    saveas(gcf,'vc_velocity.jpg');
end


function [] = Result_plots( x,y,imax,jmax,dx,v,p,t,Re)
    streamFunction=zeros(imax,jmax);
    % Streamfunction from v=d(w)/dx 
    for i=2:imax-1
       streamFunction(i,2:jmax-1)=streamFunction(i-1,2:jmax-1)-dx*v(i,2:jmax-1);
    end
    
    figure
    contourf(x,y,streamFunction',10);
    colorbar;
    axis([0 1 0 1]);
    ylabel('Y=y/L');
    xlabel('X=x/L');
    figTitle=sprintf('Stream function at t=%ss for Re=%s',num2str(t),num2str(Re));
    title(figTitle);
    figName=sprintf('SFat%ssat%sasRe.jpg',num2str(t),num2str(Re));
    saveas(gcf,figName);
   
    x_center=zeros(length(x)-1,1);
    y_center=zeros(length(y)-1,1);
    for i=1:length(x_center)
        x_center(i)=(x(i)+x(i+1))/2;
    end
    for j=1:length(y_center)
        y_center(j)=(y(j)+y(j+1))/2;
    end


    figure
    contourf(x_center,y_center,p(2:end-1,2:end-1)',10);
    colorbar;
    axis([0 1 0 1]);
    xlabel('X=x/L');
    ylabel('Y=y/L');
    figTitle=sprintf('P contours at t=%ss for Re=%s',num2str(t),num2str(Re));
    title(figTitle);
    figName=sprintf('pContourat%ssat%sasRe.jpg',num2str(t),num2str(Re));
    saveas(gcf,figName);
end

