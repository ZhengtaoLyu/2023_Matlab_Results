clear
close all
clc 
%% prey agents parameters 1
N = 20; 
T = 500; % simulation time
Re = 4; %猎物半径
Rr = 0.8; % repulsion radius
Ro = 5; % orientation radius
Ra = 10; % attraction radius
base_speed = 0.15; % speed
max_theta = 1.75; % maximum turning angle
sigma = 0.02; % standard deviation for directional noise
death_pray_distance = 0.5;

%% predator agents parameters
pN = 3; 
Rs = 6; % predator radius
pbase_speed = 0.16; % speed
pmax_theta = 1.5; % maximum turning angle
psigma = 0.01; % standard deviation for directional noise
death_predator_times = 8;

%% prey agents initialisation
C = NaN(N, T+1, 2); %C is the position matrix, the scale is 100*501*2, each agent has 2 position components, C(:,:,1) is the X coordinate, C( :,:,2) is the Y coordinate
V = NaN(N, T+1, 2); %V is the velocity matrix, the scale is 100*501*2, each agent has 2 velocity components, V(:,:,1) is the X velocity, V( :,:,2) is the Y speed
Speeds = NaN(N, T+1); %Speeds is a speed scalar, indicating the speed of the agent
Angles = NaN(N, T+1); %Angles is an angle scalar, indicating the direction of the agent
C(:,1,1) = sqrt(N)*randn(N,1,1);% x coordinate at time %1
C(:,1,2) = sqrt(N)*randn(N,1,1);% y coordinate at time %1
Angles(:, 1) = random('Uniform', -pi, pi, N, 1); %The angle at the initial moment, the angle of each agent is random
Speeds(:, 1) = base_speed*ones(N,1); %The speed at the initial moment, the speed of each agent is the same
V(:,1,:) = base_speed*[cos(Angles(:,1)) sin(Angles(:,1))]; %The velocity matrix at the initial moment, related to the initial angle

%% predator agents initialisation
pC = NaN(pN, T+1, 2); %C is the position matrix, the scale is 100*501*2, each agent has 2 position components, C(:,:,1) is the X coordinate, C( :,:,2) is the Y coordinate
pV = NaN(pN, T+1, 2); %V is the velocity matrix, the scale is 100*501*2, each agent has 2 velocity components, V(:,:,1) is the X velocity, V( :,:,2) is the Y speed
pSpeeds = NaN(pN, T+1); %Speeds is a speed scalar, indicating the speed of the agent
pAngles = NaN(pN, T+1); %Angles is an angle scalar, indicating the direction of the agent
pC(:,1,1) = sqrt(N)*randn(pN,1,1); %x coordinate at time %1
pC(:,1,2) = sqrt(N)*randn(pN,1,1); %y coordinate at time %1
pAngles(:, 1) = random('Uniform', -pi, pi, pN, 1); %The angle at the initial moment, the angle of each agent is random
pSpeeds(:, 1) = pbase_speed*ones(pN,1); %The speed at the initial moment, the speed of each agent is the same
pV(:,1,:) = pbase_speed*[cos(pAngles(:,1)) sin(pAngles(:,1))]; %The velocity matrix at the initial moment, related to the initial angle
%% simulation
for t = 1:T
    %predator and prey distance
    X = repmat(squeeze(C(:,t,1)),1,pN);
    Y = repmat(squeeze(C(:,t,2)),1,pN); 
    pX = repmat(squeeze(pC(:,t,1)),1,N); 
    pY = repmat(squeeze(pC(:,t,2)),1,N);
    ppDistance = sqrt((transpose(pX) - X).^2 + (transpose(pY) - Y).^2);

    choice = transpose(ppDistance)<Rs;
choicesum = sum(choice,2); % is 0, which means there is no prey in the field of vision, so walk at the original speed
    choicesum = choicesum>0;
    choice = choice.*transpose(ppDistance);

    ind = find(choice==0);
    choice(ind) = inf;
    [m,n] = min(choice'); %n is the location
    
    ppT = ppDistance';
    ppX = X - transpose(pX);
    ppY = Y - transpose(pY);
    pds1(:,t) = [ ppX(n(1),1)/ppT(1,n(1)) ; ppX(n(2),2)/ppT(2,n(2)) ; ppX(n(3),3)/ppT(3,n(3)) ];
    pds2(:,t) = [ ppY(n(1),1)/ppT(1,n(1)) ; ppY(n(2),2)/ppT(2,n(2)) ; ppX(n(3),3)/ppT(3,n(3)) ];
    pds = cat(3,pds1,pds2);

    pdi = ~choicesum.*pV(:,t,:) + choicesum.* pds(:,t,:);

    danger = ppDistance<Re; %100*3, 1 means a predator appears and you have to run
    dangersum =  sum(danger, 2);
    danger_distances = danger.*ppDistance; 
    danger_X_pos_differences = danger.*(transpose(pX) - X); 
    danger_Y_pos_differences = danger.*(transpose(pY) - Y); 
    de = -(sum(cat(3,danger_X_pos_differences./danger_distances,danger_Y_pos_differences./danger_distances), 2, "omitnan"));
    
    %prey and prey distance
    X = repmat(squeeze(C(:,t,1)),1,N); %squeeze(C(:,t,1)) is 100*1, X coordinate.
    Y = repmat(squeeze(C(:,t,2)),1,N);%squeeze(C(:,t,2)) is 100*1, Y coordinate.
    distances = sqrt((transpose(X) - X).^2 + (transpose(Y) - Y).^2); %Symmetric matrix, ij is the distance between agent i and agent j
    distances = distances + diag(Inf.*ones(1,N)); %Symmetric matrix, diag(Inf.*ones(1,N)) is a diagonal matrix whose diagonal is Inf
    
    %
    repels = distances < Rr;% Symmetric matrix, the part of 1 means that the ij distance is smaller than Rr, and the diagonal is all 0
    orients = (distances < Ro).*(distances >= Rr); % Symmetric matrix, the part of 1 means that the ij distance is greater than Rr and less than Ro, and the diagonal is all 0
    attracts = (distances < Ra).*(distances >= Ro); % Symmetric matrix, the part of 1 means that the ij distance is greater than Ro and less than Ra, and the diagonal is all 0

    X_pos_differences = transpose(X) - X; %Antisymmetric matrix, about the X coordinate
    Y_pos_differences = transpose(Y) - Y; %Antisymmetric matrix, about the Y coordinate
    repel_distances = repels.*distances; %The diagonal is all 0
    repel_X_pos_differences = repels.*X_pos_differences;
    repel_Y_pos_differences = repels.*Y_pos_differences;
    dr = -(sum(cat(3,repel_X_pos_differences./repel_distances,repel_Y_pos_differences./repel_distances), 2, "omitnan"));%'omitnan')Calculate the elemental sum of NaN values, dr is 100*1*2
    dr = ~dangersum.*dr;
    %
    orientations = squeeze(V(:,t,:)); %orientations are all velocity vectors
    neighbour_speeds = sqrt(orientations(:,1).^2 + orientations(:,2).^2); %neighbour_speedsare all velocity scalars
    norm_orientations = orientations./repmat(neighbour_speeds, 1, 2); %100*2
    do = sum(repmat(reshape(norm_orientations, 1, N, 2), N, 1, 1).*orients, 2, "omitnan");
    %
    attract_distances = attracts.*distances;
    attract_X_pos_differences = attracts.*X_pos_differences;
    attract_Y_pos_differences = attracts.*Y_pos_differences;
    da = (sum(cat(3, attract_X_pos_differences./attract_distances,attract_Y_pos_differences./attract_distances), 2, "omitnan"));
    %
    no_repel = sum(repels,2) == false; %The all-zero row in %repels, the i-th row in no_repel is 1, which means agent i is far away from other agents
    no_repel = ~dangersum.*no_repel;
    no_interact = sum(repels,2) + sum(orients,2) + sum(attracts,2) == false; %The i-th line in %no_interact is 1, indicating that there are no other agents within the Ra range around agent i.
    no_interact = ~dangersum.*no_interact;
    
    di = de + dr + (0.5*do + 0.5*da).*no_repel + V(:,t,:).*no_interact;

    di_angle = angle(di(:,:,1) + 1i*di(:,:,2)) + random('Normal', 0, sigma, N, 1);
    theta = wrapToPi(di_angle - Angles(:,t));
    theta = min(theta, max_theta);
    theta = max(theta, -max_theta);
    %
    Speeds(:,t+1) = Speeds(:,t);
    Angles(:,t+1) = Angles(:,t) + theta;
    V(:,t+1,:) = Speeds(:,t+1).*[cos(Angles(:,t+1)), sin(Angles(:,t+1))];
    C(:,t+1,:) = C(:,t,:) + V(:,t+1,:);

    pdi_angle = angle(pdi(:,:,1) + 1i*pdi(:,:,2)) + random('Normal', 0, psigma, pN, 1);
    ptheta = wrapToPi(pdi_angle - pAngles(:,t));
    ptheta = min(ptheta, pmax_theta);
    ptheta = max(ptheta, -pmax_theta);
    %
    pSpeeds(:,t+1) = pSpeeds(:,t);
    pAngles(:,t+1) = pAngles(:,t) + ptheta;
    pV(:,t+1,:) = pSpeeds(:,t+1).*[cos(pAngles(:,t+1)), sin(pAngles(:,t+1))];
    pC(:,t+1,:) = pC(:,t,:) + pV(:,t+1,:);
end
%% visualisation/animation
minx = min(C(:,:,1), [], "all");
maxx = max(C(:,:,1), [], "all");
miny = min(C(:,:,2), [], "all");
maxy = max(C(:,:,2), [], "all");
%
for t = 0:T
    quiver(C(:,t+1,1), C(:,t+1,2), V(:,t+1,1), V(:,t+1,2) ,'b')
    hold on
    quiver(pC(:,t+1,1), pC(:,t+1,2), pV(:,t+1,1), pV(:,t+1,2) ,'r')
    
    axis equal
    axis([minx, maxx, miny, maxy])
    title("t = " + t)
    hold off
    drawnow
    pause(0.5)
end