function simulation()

function surface()
    lighting phong;
    material shiny;
    lightangle(-45,30)
    light('Position',[-10 20 10]);
    axis([1 60 1 60 -2 8]);
end


function newH=Wave(n,theta,rho,dt,c,k,H,oldH,fix,cont,connect)
    global Ekin Epot;
         
    % DAMPED WAVE EQUATION:
    %
    % d^2/dt^2*h + K*(dh/dt) = C^2*(d^2*h/dx^2 + d^2*h/dy^2)
    %
    %   where   h := Height
    %           K := Damping Constant
    %           C := Wave Speed
    %   The right side of the equation is the potential (height of one
    %   element regarding its neighbours).
    %   The wave equation implies that acceleration (d^2*h/dt^2) and 
    %   velocity (dh/dt) of each element are produced through its potential.
    %
    % Finite Difference Procedure:
    %
    %   velocity     := dh/dt      :=   (H(i,j)-oldH(i,j))/dt
    %   acceleration := d^2/dt^2*h :=   ((newH(i,j)-H(i,j))-(H(i,j)-oldH(i,j)))/dt^2
    %                               =   (newH(i,j)-2*H(i,j)+oldH(i,j))/dt^2
    %   similar we have:
    %
    %   potential in x-direction := d^2/dx^2*h := (H(i+1,j)-2*H(i,j)+H(i-1,j))/dx^2
    %   potential in y-direction := d^2/dy^2*h := (H(i,j+1)-2*H(i,j)+H(i,j-1))/dy^2
    %   where dx=dy=1; (spacing between 2 points);
    %   It is possible to include the potential in DIAGONAL direction as well.
    %   Nummerical experiments demomstrate that wave shapes look visually
    %   better when all 8 neighbours are taken into account.
    %   But remember, dx and dy, (spacing between 2 points) in diagonal direction is
    %   sqrt(2) and not 1. Therefor 1/dx^2 resp. 1/dy^2 is equal 0.5
    
    %   Apply these to the wave equation above we have:
    
    potential = zeros(length(theta),length(rho));
    velocity = zeros(length(theta),length(rho));
    acceleration = zeros(length(theta),length(rho));
    newH = zeros(length(theta),length(rho));
    for t = 1:length(theta)
        for r = 1:length(rho)
            if t == 1 && r == 1
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(mod(length(theta)/2+t-1,length(theta))+1,r))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(t+1,r)-2*H(t,r)+H(length(theta),r))));
            elseif t == length(theta) && r == 1
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(mod(length(theta)/2+t-1,length(theta))+1,r))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(1,r)-2*H(t,r)+H(t-1,r))));
            
            elseif t == 1 && r == length(rho)
                potential(t,r) = 0;
            elseif t == length(theta) && r == length(rho)
                potential(t,r) = 0;
                
            
            elseif t == 1
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(t,r-1))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(t+1,r)-2*H(t,r)+H(length(theta),r))));
            elseif t == length(theta)
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(t,r-1))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(1,r)-2*H(t,r)+H(t-1,r))));
                
                
            elseif r == 1
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(mod(length(theta)/2+t-1,length(theta))+1,r))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(t+1,r)-2*H(t,r)+H(t-1,r))));
                
            elseif r == length(rho)
                potential(t,r)=0;
                
            else
                potential(t,r) = c^2*((H(t,r+1)-2*H(t,r)+H(t,r-1))+((1./r)*(H(t,r+1)-H(t,r)))+((1./r.^2)*(H(t+1,r)-2*H(t,r)+H(t-1,r))));
            end
            velocity(t,r)=(H(t,r)-oldH(t,r))/dt;
            acceleration(t,r)=-k*velocity(t,r)+potential(t,r);
            newH(t,r)=acceleration(t,r)*dt^2-oldH(t,r)+2*H(t,r);
        end
    end
    
    
%     -c^2*((4*H(theta,rho)-H(theta+1,rho)-H(theta-1,rho)-H(theta,rho+1)-H(theta,rho-1))...
%         +0.5*(4*H(theta,rho)-H(theta+1,rho+1)-H(theta+1,rho-1)-H(theta-1,rho+1)-H(theta-1,rho-1))); %  diagonal direction (opitonal)
%     velocity(theta,rho)=(H(theta,rho)-oldH(theta,rho))/dt;
%     acceleration(theta,rho)=-k*velocity(theta,rho)+potential(theta,rho); %  := (newH(i,j)-2*H(i,j)+oldH(i,j))/dt^2 as mentioned above
    % therefor, the new height is:
%     newH(theta,rho)=acceleration(theta,rho)*dt^2-oldH(theta,rho)+2*H(theta,rho);
      
    % Please take notice that this equation isn't applied for the elements
    % along the edges and at the corners (Boundary Points / Randpunkte),
    % that's why i and j are from 2 to n-1 instead of 1 to n.
    
    
%     
%     % BOUNDARY CONDITIONS:
%     %
%     %   Equations for boundary points.
%     %   Keep in mind that elements along the edges have 5 neighbours
%     %   instead of 8 and vertices only have 3.
%     
%     
%     
%     potential(n,rho)=-c^2*((3*H(n,rho)-H(n-1,rho)-H(n,rho+1)-H(n,rho-1))...
%         +0.5*(2*H(n,rho)-H(n-1,rho+1)-H(n-1,rho-1)));
%     potential(theta,n)=-c^2*((3*H(theta,n)-H(theta,n-1)-H(theta+1,n)-H(theta-1,n))...
%         +0.5*(2*H(theta,n)-H(theta+1,n-1)-H(theta-1,n-1)));
%     potential(1,rho)=-c^2*((3*H(1,rho)-H(2,rho)-H(1,rho+1)-H(1,rho-1))...
%         +0.5*(2*H(1,rho)-H(2,rho+1)-H(2,rho-1)));
%     potential(theta,1)=-c^2*((3*H(theta,1)-H(theta,2)-H(theta+1,1)-H(theta-1,1))...
%         +0.5*(2*H(theta,1)-H(theta+1,2)-H(theta-1,2)));
%     velocity(n,rho)=(H(n,rho)-oldH(n,rho))/dt;
%     velocity(theta,n)=(H(theta,n)-oldH(theta,n))/dt;
%     velocity(1,rho)=(H(1,rho)-oldH(1,rho))/dt;
%     velocity(theta,1)=(H(theta,1)-oldH(theta,1))/dt;
%     
%                 % 4 corners:
%     potential(1,1)=-c^2*((2*H(1,1)-H(2,1)-H(1,2))...
%         +0.5*(H(1,1)-H(2,2)));
%     potential(1,n)=-c^2*((2*H(1,n)-H(1,n-1)-H(2,n))...
%         +0.5*(H(1,n)-H(2,n-1)));
%     potential(n,1)=-c^2*((2*H(n,1)-H(n,2)-H(n-1,1))...
%         +0.5*(H(n,1)-H(n-1,2)));
%     potential(n,n)=-c^2*((2*H(n,n)-H(n-1,n)-H(n,n-1))...
%         +0.5*(H(n,n)-H(n-1,n-1)));
%     velocity(1,1)=(H(1,1)-oldH(1,1))/dt;
%     velocity(1,n)=(H(1,n)-oldH(1,n))/dt;
%     velocity(n,1)=(H(n,1)-oldH(n,1))/dt;
%     velocity(n,n)=(H(n,n)-oldH(n,n))/dt;
% 
% 
% 
% %                       DEFAULT MODUS     
% if(~fix && ~cont && ~connect)
%     
%     acceleration(n,rho)=-k*velocity(n,rho) +potential(n,rho);
%     newH(n,rho)=acceleration(n,rho)*dt^2-oldH(n,rho)+2*H(n,rho);
%     
%     acceleration(theta,n)=-k*velocity(theta,n) +potential(theta,n);
%     newH(theta,n)=acceleration(theta,n)*dt^2-oldH(theta,n)+2*H(theta,n);
%     
%     acceleration(1,rho)=-k*velocity(1,rho) +potential(1,rho);
%     newH(1,rho)=acceleration(1,rho)*dt^2-oldH(1,rho)+2*H(1,rho);
%     
%     acceleration(theta,1)=-k*velocity(theta,1) +potential(theta,1);
%     newH(theta,1)=acceleration(theta,1)*dt^2-oldH(theta,1)+2*H(theta,1);
%     
%     
%                 % 4 corners:
%     acceleration(1,1)=-k*velocity(1,1) +potential(1,1);
%     newH(1,1)=acceleration(1,1)*dt^2-oldH(1,1)+2*H(1,1);
%     
%     acceleration(1,n)=-k*velocity(1,n) +potential(1,n);
%     newH(1,n)=acceleration(1,n)*dt^2-oldH(1,n)+2*H(1,n);
%     
%     acceleration(n,1)=-k*velocity(n,1) +potential(n,1);
%     newH(n,1)=acceleration(n,1)*dt^2-oldH(n,1)+2*H(n,1);
%     
%     acceleration(n,n)=-k*velocity(n,n) +potential(n,n);
%     newH(n,n)=acceleration(n,n)*dt^2-oldH(n,n)+2*H(n,n);   
%     
% end

    kin=velocity.^2;
    pot=-acceleration.*oldH;
    Ekin=[Ekin sum(kin(:))];
    Epot=[Epot sum(pot(:))];
end

    n = 60;
    H=zeros(n,n);
    H(1,:) = 1;
    H(2,:) = 1;
    
    oldH=H;
    newH=H;
    theta = 0:2*pi/(n-1):2*pi; 
    rho = 0:n-1;
    [r,t] = meshgrid(rho,theta);
    [x,y,v] = pol2cart(t,r,newH);
    figure('color','white');
    [Xi,Yi,Zi,g] = polarplot3d(v);

%     h=surf(x,y,v(2:n-1,2:n-1));
%     surface();
    

    while true
%         run=get(handles.RUN_button,'value');
%         if ~run
%             break
%         end
        newH=Wave(n,theta,rho,0.05,1,0,H,oldH,0,0,0);
                % n,i,j,dt,c,k,H,oldH,fix,cont,connect
        
        [x,y,v] = pol2cart(t,r,newH);
                
        set(g,'zdata',v);
        pause(0.05);
        oldH=H;
        H=newH;
    end


end


