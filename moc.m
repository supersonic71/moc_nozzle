%A method of characteristics implementation for min-length nozzle

%All angles stored in degrees
%flowprandtlmeyer i/o is degree
%slope_inter expects degrees

m_exit = input("Enter your exit Mach number : ");
figure;
gamma = 1.4;
[~, nu_exit, mu] = flowprandtlmeyer(gamma, m_exit, 'mach');
th_max = nu_exit./2;

fan_start = .3;
n_waves = 70;
fan_th = linspace(fan_start,th_max,n_waves); 

n_points = (n_waves+2).*(n_waves+1)./2 - 1; %Sum of first (1+n_waves) natural numbers minus one



%symbols
KM = 1;
KP = 2;
TH = 3;
NU = 4;
MM = 5;
MU = 6;
X=7;
Y=8;
n_properties = 8;
grid = zeros(n_points, n_properties);

%checking error
%grid(n_points,:) = [th_max*2, -th_max*2, 0, th_max*2, m_exit,mu];

%can we assume tan(th) = th
r = 2.5; %throat radius
%first hits of expansion waves
for n_point = 1:n_waves
    th = fan_th(n_point);
    nu = th;
    Km = th+nu;
    Kp = th-nu;
    [M, ~, mu] = flowprandtlmeyer(gamma, nu, 'nu');
    
    if n_point == 1
        m1 = 0;
        x1 = 0;
        y1 = 0;
    else
        m1 = mean([th,grid(n_point-1,TH)]) + mean([mu,grid(n_point-1,MU)]); %left running characteristic
        x1=grid(n_point-1,X);
        y1=grid(n_point-1,Y);
    end
    
    m2 = th-mu;
    x2=0;
    y2=r;
    [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
    
    
    grid(n_point,:) = [Km, Kp, th, nu, M, mu, xo, yo];
end

%test
%plot(grid(1:7,[7])' , grid(1:7,[8])')


%th = 0 condition
%then rest are intersection condition
%copy previous to the one hitting wall
%repeat
n_point = n_point + 1;
grid(n_point,:) = grid(n_point-1,:);
x1=grid(n_point-1,X);
y1=grid(n_point-1,Y);
m1=grid(n_point,TH) + grid(n_point,MU); 
x2=0;
y2=r;
m2=mean([grid(n_point,TH),th_max]); 
[xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
%fprintf("%f ",x1,y1,m1,x2,y2,m2,xo,yo); 
grid(n_point,[X,Y]) = [xo, yo];

for n_wave = 2:n_waves
    
    
    
    %th = 0 condition
    n_point = n_point+1;
    Km = grid(n_wave, KM);
    Kp = -Km;
    th = 0;
    nu = .5.*(Km-Kp); 
    [M, ~, mu] = flowprandtlmeyer(gamma, nu, 'nu');
    grid(n_point,:) = [Km, Kp, th, nu, M, mu,0,0];
    
    m1 = 0;
    x1 = 0;
    y1 = 0;
    cn = n_waves-n_wave+2;
    x2=grid(n_point-cn,X);
    y2=grid(n_point-cn,Y);
    m2=mean([grid(n_point,TH),grid(n_point-cn,TH)]) - mean([grid(n_point,MU),grid(n_point-cn,MU)]);
    [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
    grid(n_point,[X,Y]) = [xo, yo];
    
    %intersections
    for i = n_wave+1:n_waves
        n_point = n_point+1;
        Km = grid(i,KM);
        Kp = grid(n_point-1, KP);
        th = .5.*(Km+Kp);
        nu = .5.*(Km-Kp);
        [M, ~, mu] = flowprandtlmeyer(gamma, nu, 'nu');
        grid(n_point,:) = [Km, Kp, th, nu, M, mu,0,0];
        
        x1=grid(n_point-1,X);
        y1=grid(n_point-1,Y);
        m1=mean([grid(n_point,TH),grid(n_point-1,TH)]) + mean([grid(n_point,MU),grid(n_point-1,MU)]); 
        cn = n_waves-n_wave+2;
        x2=grid(n_point-cn,X);
        y2=grid(n_point-cn,Y);
        m2=mean([grid(n_point,TH),grid(n_point-cn,TH)]) - mean([grid(n_point,MU),grid(n_point-cn,MU)]);
        [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
        grid(n_point,[X,Y]) = [xo, yo];
    
    end
    
    
    cn = n_waves-n_wave+2;
    %copy
    n_point = n_point + 1;
    grid(n_point,:) = grid(n_point-1,:);
    
    x1=grid(n_point-1,X);
    y1=grid(n_point-1,Y);
    m1=grid(n_point,TH) + grid(n_point,MU); 
    x2=grid(n_point-cn,X);
    y2=grid(n_point-cn,Y);
    m2=mean([grid(n_point,TH), grid(n_point-cn,TH)]); 
    [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
    grid(n_point,[X,Y]) = [xo, yo];
end
% 
% plot(grid(:,[7])' , grid(:,[8])','r');
% hold on;
% plot([0 grid(n_waves+1,X)],[r grid(n_waves+1,Y)],'r');
% plot(grid(:,[7])' , -grid(:,[8])','r');
% hold on;
% plot([0 grid(n_waves+1,X)],-[r grid(n_waves+1,Y)],'r');
% axis equal
contour_points = zeros(1,n_waves);
contour_points(1) = n_waves + 1;
for n_wave = 2:n_waves
    contour_points(n_wave) = contour_points(n_wave-1) + (n_waves+2-n_wave);
end
hold on;
plot([0 grid(contour_points,X)'],[r grid(contour_points,Y)'],'r','LineWidth',3);
plot([0 grid(contour_points,X)'],-[r grid(contour_points,Y)'],'r','LineWidth',3);
axis equal
title("Best minimum length nozzle for Mach " + m_exit);

