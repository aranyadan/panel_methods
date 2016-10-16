clear all; 
clc;
Qinf = 10;
hold on;

%% Vortex panel method with a constant source sheet

%  Open a File and read airfoil coordinates
%
fid = fopen('naca_0015.dat','r');
n = fscanf(fid,'%d',1);



x = zeros(n);
y = zeros(n);

for i=1:n
    a=fscanf(fid,'%f',1);
    b=fscanf(fid,'%f',1);
    x(n-i+1) = a;
    y(n-i+1) = b;
end
fclose(fid);
alpha = 0;
%
%% Plot the airfoil on window #1
%
% plot(x,y);
n=length(x) - 1; %%n is the no of panels
A=zeros(n,n);
ds=zeros(1,n);
alphai=zeros(1,n);
collocpoints = zeros(2,n);  %% 1 is x, 2 is y



for i=1:n
    alphai(i) = atan2d(y(i) - y(i+1) , x(i) -  x(i+1));
    if (alphai(i) > 0)
      alphai(i) = 180 - alphai(i);
    else
      alphai(i) = -180 - alphai(i);
    end    
    ds(i) = sqrt((x(i) - x(i+1))^2 + (y(i) - y(i+1))^2 );
    collocpoints(1,i) = (x(i) + x(i+1))/2;
    collocpoints(2,i) = (y(i) + y(i+1))/2;   
end
%% Influence coefficient matrix A
for i=1:n
    for j=1:n
        if i==n
            A(i,1) = 1;
            A(i,n) = 1;
            i=i+1;
            break;
        end
        if i==j
            A(i,j) = -0.5;
            continue;
        end
        [up,vp] = VOR2DC(1,collocpoints(1,i),collocpoints(2,i),x(j),y(j),x(j+1),y(j+1));
        Rot = [cosd(alphai(j)) sind(alphai(j)) ; -sind(alphai(j)) cosd(alphai(j))];
        pVEL = [up;vp];
        VEL = Rot * pVEL;
        u = VEL(1);
        v = VEL(2);

        A(i,j) = u * cosd(alphai(i)) - v * sind(alphai(i));
    end
end
%% For source part

B=A;
B(:,n+1)=0;
for i=1:n
    for j=1:n
        if i==j
            continue;
        end
        [up,vp] = SOR2DC(1,collocpoints(1,i), collocpoints(2,i),x(j),y(j),x(j+1),y(j+1));
        Rot = [cosd(alphai(j)) sind(alphai(j)) ; -sind(alphai(j)) cosd(alphai(j))];
        pVEL = [up;vp];
        VEL = Rot * pVEL;
        u = VEL(1);
        v = VEL(2);

        B(i,n+1) = B(i,n+1) + u * cosd(alphai(i)) - v * sind(alphai(i));
    end
end
B(n+1,1) = 1;
B(n+1,n) = 1;
alpha = 1:10;
%% Setting up RHS
for k = 1:10
    RHS = zeros(n+1,1);
    for i=1:n
        Uinf = Qinf * cosd(alpha(k));
        Vinf = Qinf * sind(alpha(k));
        RHS(i,1) = Vinf * sind(alphai(i)) - Uinf * cosd(alphai(i));
    end
    RHS(n+1,1) = 0;
    gamma = RHS\B;

    Cp = zeros(1,n);
    for i=1:n
       Cp(i) = 1 -  ((Qinf * cosd(alpha(k) + alphai(i)) + gamma(i)/2)/Qinf)^2;
        
    end
    for i=1:n
        dl(i) =2*gamma(i) * ds(i) / Qinf ;
    end
    L(k)=sum(dl);
end



plot(collocpoints(1,:),-Cp);
% plot(alpha,L);