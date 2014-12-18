function simulation()

function surface()
    lighting phong;
    material shiny;
    lightangle(-90,90)
    light('Position',[-1 2 1]);
    axis([-1.2 1.2 -1.2 1.2 -0.5 0.5]);
    grid off;
end


function [newH,newCent]=Wave(n,rho,theta,dt,c,k,H,oldH,Cent,fix,cont,connect)
    global Ekin Epot;
    
    % The original cartesian definition of this function we repurposed to
    % calculate the polar discretized model.
         
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
    
    potential = zeros(length(rho),length(theta));
    velocity = zeros(length(rho),length(theta));
    acceleration = zeros(length(rho),length(theta));
    newH = zeros(length(rho),length(theta));
    for t = 1:length(theta)
        for r = 1:length(rho)
            if t == 1 && r == 1
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+Cent)+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,t+1)-2*H(r,t)+H(r,length(theta)))/(2*pi/length(theta))^2));
            elseif t == length(theta) && r == 1
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+Cent)+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,1)-2*H(r,t)+H(r,t-1))/(2*pi/length(theta))^2));
            
            elseif t == 1 && r == length(rho)
                potential(r,t) = 0;
            elseif t == length(theta) && r == length(rho)
                potential(r,t) = 0;
                
            
            elseif t == 1
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+H(r-1,t))+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,t+1)-2*H(r,t)+H(r,length(theta)))/(2*pi/length(theta))^2));
            elseif t == length(theta)
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+H(r-1,t))+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,1)-2*H(r,t)+H(r,t-1))/(2*pi/length(theta))^2));
                
                
            elseif r == 1
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+Cent)+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,t+1)-2*H(r,t)+H(r,t-1))/(2*pi/length(theta))^2));
                
            elseif r == length(rho)
                potential(r,t)=0;
                
            else
                potential(r,t) = c^2*((H(r+1,t)-2*H(r,t)+H(r-1,t))+((1./r)*(H(r+1,t)-H(r,t)))+((1./r.^2)*(H(r,t+1)-2*H(r,t)+H(r,t-1))/(2*pi/length(theta))^2));
            end
            velocity(r,t)=(H(r,t)-oldH(r,t))/dt;
            acceleration(r,t)=-k*velocity(r,t)+potential(r,t);
            newH(r,t)=acceleration(r,t)*dt^2-oldH(r,t)+2*H(r,t);
        end
    end
    newCent=mean(newH(1,:));

    kin=velocity.^2;
    pot=-acceleration.*oldH;
    Ekin=[Ekin sum(kin(:))];
    Epot=[Epot sum(pot(:))];
end

    n = 60;
    m = 60;
    H=zeros(n,m);
    
    % Radially Symmetric Gaussian Initial Condition
%     sigma = 3;    % just an example value
%     gSize = 3*sigma;  % cutoff point
%     xNum = 0:gSize;
%     H(1,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-1.^2 / (2*sigma^2)));
%     H(2,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-2.^2 / (2*sigma^2)));
%     H(3,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-3.^2 / (2*sigma^2)));
%     H(4,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-4.^2 / (2*sigma^2)));
%     H(5,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-5.^2 / (2*sigma^2)));
%     H(6,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-6.^2 / (2*sigma^2)));
%     H(7,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-7.^2 / (2*sigma^2)));
%     H(8,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-8.^2 / (2*sigma^2)));
%     H(9,:) = 6*(1 / (sigma * sqrt(2 * pi)) * exp(-9.^2 / (2*sigma^2)));

    % A (mediocre) attempt at a non-centered smooth peak
    H(25,30) = .5;
    H(26,30) = .45;
    H(27,30) = .35;
    H(28,30) = .25;
    H(29,30) = .15;
    H(30,30) = .1;
    H(31,30) = .075;
    H(32,30) = .06;
    H(33,30) = .05;
    H(34,30) = .04;
    H(35,30) = .03;
    H(36,30) = .02;
    H(37,30) = .01;
    
    H(24,30) = .45;
    H(23,30) = .35;
    H(22,30) = .25;
    H(21,30) = .15;
    H(20,30) = .1;
    H(19,30) = .075;
    H(18,30) = .06;
    H(17,30) = .05;
    H(16,30) = .04;
    H(15,30) = .03;
    H(14,30) = .02;
    H(13,30) = .01;
    
    H(25,31) = .3;
    H(25,32) = .1;
    H(25,33) = .05;
    H(25,34) = .01;
    
    H(25,29) = .3;
    H(25,28) = .1;
    H(25,27) = .05;
    H(25,26) = .01;
    
    H(24,31) = .25;
    H(24,32) = .125;
    
    H(23,31) = .125;
    H(23,32) = .02;
    
    H(24,29) = .25;
    H(24,28) = .125;
    
    H(23,29) = .125;
    H(23,28) = .02;
    
    H(26,31) = .25;
    H(26,32) = .125;
    
    H(27,31) = .125;
    H(27,32) = .02;

    H(26,29) = .25;
    H(26,28) = .125;
    
    H(27,29) = .125;
    H(27,28) = 0.02;

    Cent = 0;
    newCent=Cent;
    
    oldH=H;
    newH=H;
    theta = 0:2*pi/m:2*pi;
    theta = theta(1:end-1);
    rho = 0:n-1;
    [rh,th] = meshgrid(rho,theta);
    [x,y,v] = pol2cart(th,rh,newH);
    figure('color','white');
    [Xi,Yi,Zi,g] = polarplot3d(v);
    surface;
    while true
        [newH,newCent]=Wave(n,rho,theta,0.05,15,.2,H,oldH,Cent,0,0,0);
                % n,i,j,dt,c,k,H,oldH,fix,cont,connect
        
        [x,y,v] = pol2cart(th,rh,newH);
        [Xi,Yi,Zi,g] = polarplot3d(v);
        surface();
        pause(0.05);
        oldH=H;
        H=newH;
        Cent = newCent;
    end


end


