clear all; 
clc;
Qinf = 10;
hold on;

%% Input
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

%% Plot the airfoil on window #1
%
plot(x,y); 
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
for i=1:n-1
  for j = 1:n
    if(i == j)
      A(i,i) = -0.5;
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
A(n,1) = 1;
A(n,n) = 1;
%% RHS:
RHS = zeros(n,1);

for i = 1:n-1
    Uinf = Qinf * cosd(alpha);
    Vinf = Qinf * sind(alpha);
    RHS(i,1) = Vinf * sind(alphai(i)) - Uinf * cosd(alphai(i));
end
gamma = A\RHS;

Cp = zeros(1,n);
for i=1:n-1
   Cp(i) = 1 -  ((Qinf * cosd(alpha + alphai(i)) + gamma(i)/2)/Qinf)^2;
   if(Cp(i) < -30)
       Cp(i) = 0;
   end
end

plot(collocpoints(1,:),Cp);

