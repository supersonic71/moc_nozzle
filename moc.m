% Ephraim Raj, K.R.S Nandhan
% A method of characteristics implementaion to find
% the proper contour for a shock free expansion of
% minimum length nozzle - rocket type as opposed to wind tunnel type

% User inputs - mach number, gamma (specific heat ratio),
%               number of expansion waves
% Output -      a file with the required coordinates of contour (TODO)
%               a figure displaying the contour
           
% All angles stored in degrees
% flowprandtlmeyer i/o is degree
% slope_inter expects degrees

% NOTE : the function flowprandtlmeyer() requires MATLAB Aerospace Toolbox
% NOTE : the slope_inter() function is required
% Report issues/feedback : github.com/sgt-miller/moc_nozzle

% Some formulas were deduced using the figure attached with Example 11.1
% of Modern Compressible Flow, 2nd Edition, Anderson Jr. 
% Tables, Figures, Equations referenced are from the above book
% Throat corner is point a in Figure 11.12

m_exit = input("Enter your exit Mach number : ");
gamma = 1.4;
[~, nu_exit, mu] = flowprandtlmeyer(gamma, m_exit, 'mach');
th_max = nu_exit./2; % Proof in Ch 11, Modern Compressible flow, Anderson Jr. 

% Generating expansion waves from throat corner
fan_start = .375; % Initial angle
n_waves = 7; % Number of expansion waves to be generated
fan_th = linspace(fan_start,th_max,n_waves); 

n_points = (n_waves+2).*(n_waves+1)./2 - 1; % Sum of first (1+n_waves) natural numbers minus one
                                            % deduced from analysing example
                                            
% For example, for 70 n_waves, we get 2555 n_points


% The matrix grid stores all properties related to all the points
% After calculation, to access a point 'n' and its property, 
% use grid(n,PROPERTY_SYMBOL)
% For example, to get the Mach angle at point 9, use grid(9,MU)

KM = 1; % Negative (right running) characteristic line
KP = 2; % Positive (left running) characteristic line
TH = 3; % Theta - angle of velocity vector wrt horizontal
NU = 4; % The Prandtl Meyer function
MM = 5; % Mach Number
MU = 6; % Mach Angle
X=7;    % X Coordinate
Y=8;    % Y Coordinate
n_properties = 8;
grid = zeros(n_points, n_properties);




r = 2.5; % throat radius

% first intersections of the expansion waves
% A new point is formed at the intersection of two lines
% To get location of new point, you need two previous points and slopes of those
% lines

for n_point = 1:n_waves 
    th = fan_th(n_point); % Refer Example 11.1 on why th at corner is the same as the first intersection
    nu = th;              % Section 11.7
    Km = th+nu;           
    Kp = th-nu;
    [M, ~, mu] = flowprandtlmeyer(gamma, nu, 'nu');
    
    if n_point == 1 % This intersects with the y=0 line
        m1 = 0;
        x1 = 0;
        y1 = 0;
    else            % These intersect with the previous line
        m1 = mean([th,grid(n_point-1,TH)]) + mean([mu,grid(n_point-1,MU)]); %left running characteristic
        x1=grid(n_point-1,X);
        y1=grid(n_point-1,Y);
    end
    
    % All these points come from the throat corner
    m2 = th-mu; 
    x2=0;
    y2=r;
    [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
    
    grid(n_point,:) = [Km, Kp, th, nu, M, mu, xo, yo];
    
end




% The first point of contour wall (apart from the throat corner itself)
n_point = n_point + 1;
grid(n_point,:) = grid(n_point-1,:); % "Same as point 7" in Table 11.1

x1=grid(n_point-1,X);
y1=grid(n_point-1,Y);
m1=grid(n_point,TH) + grid(n_point,MU); 

x2=0;
y2=r;
m2=mean([grid(n_point,TH),th_max]); 

[xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
grid(n_point,[X,Y]) = [xo, yo];


% From now on, we will encounter in a loop
% One point lying in the symmetric axis
% a number of points formed by intersection of expansion waves
% a point which is part of nozzle contour
% in that order
% In each loop, one expansion wave and its reflection is tracked


for n_wave = 2:n_waves
    
    
    
    % Point lying on line of symmetry
    n_point = n_point+1;
    
    Km = grid(n_wave, KM);
    th = 0; % Lies on line of symmetry
    Kp = -Km; 
    nu = .5.*(Km-Kp); 
    [M, ~, mu] = flowprandtlmeyer(gamma, nu, 'nu');
    grid(n_point,:) = [Km, Kp, th, nu, M, mu,0,0];
    
    % Finding geometric location ( y coordinate is zero, x coordinate is
    % unknown)
    
    % Intersects with y=0 line
    m1 = 0;
    x1 = 0;
    y1 = 0;
    
    
    cn = n_waves-n_wave+2; % means count of n, formula found using
                           % pattern identification in Figure 11.12
    
    x2=grid(n_point-cn,X);
    y2=grid(n_point-cn,Y);
    m2=mean([grid(n_point,TH),grid(n_point-cn,TH)]) - mean([grid(n_point,MU),grid(n_point-cn,MU)]);
    
    [xo,yo] = slope_inter([x1,x2], [y1,y2], m1, m2);
    grid(n_point,[X,Y]) = [xo, yo];
    
    % Points generated due to intersection of expansion waves
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
    
    % The point on the contour
    cn = n_waves-n_wave+2;
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


% TODO : assert the following condition
% grid(n_points,:) = [th_max*2, -th_max*2, 0, th_max*2, m_exit,mu];  


% Visualisation and output

contour_points = zeros(1,n_waves); % This has the points which are part of contour
                                   % Apart from throat corner
contour_points(1) = n_waves + 1;
for n_wave = 2:n_waves
    contour_points(n_wave) = contour_points(n_wave-1) + (n_waves+2-n_wave);
end

% Contour visualisation
figure;
hold on;
plot([0 grid(:,X)'],[r grid(:,Y)'],'r*');
plot([0 grid(:,X)'],-[r grid(:,Y)'],'r*'); %Symmetry
axis equal
title("MOC minimum length nozzle expansion section for Mach " + m_exit);

% Contour coordinates in file
% These can be exported to any modelling tools which supports importing
% coordinate files



% Note : writematrix was introduced in R2019a
% contour_coords = [grid(contour_points,X) grid(contour_points,Y)];
% writematrix(contour_coords, "contour_coords.txt", "Delimiter", "comma");

% For older versions
fid = fopen("contour_coords.txt","wt");

fprintf(fid, '%f, %f\n', 0, r); %Throat corner
for ii = contour_points
    fprintf(fid, '%f, %f\n', grid(ii,X), grid(ii,Y)); %Throat corner
end
fclose(fid);


