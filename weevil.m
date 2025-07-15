clear
close all
clc
%% part1  initial value
P = [0, 0, 0, 0, 21, 115; 
    0.32, 0.08, 0, 0, 0, 0; 
    0, 0.31, 0.11, 0, 0, 0; 
    0, 0, 0.41, 0.39, 0, 0; 
    0, 0, 0, 0.16, 0.44, 0; 
    0, 0, 0, 0, 0.21, 0.42];
eigs(P,1);

N = [6204; 2335; 823; 482; 119; 37];%initial number of weevil
%% iteration
Nn = N;
for n = 1:100
    N = P*N; 
    Nn = [Nn, N]; % save the result of every iteration
end

%% plot
plot(0:100,Nn,'linewidth',2)
xlabel('Time [days]')
ylabel('Population count [individuals]')
legend('Eggs','Neonate larvae','Later-instar larvae','Pupae','Juveniles','Adults')
   












 
%% part2 change the eigenvalue
clear
close all
clc
P = [0, 0, 0, 0, 21, 115;
0.32, 0.08, 0, 0, 0, 0;
0, 0.31, 0.11, 0, 0, 0;
0, 0, 0.41, 0.39, 0, 0;
0, 0, 0, 0.16-0.0502, 0.44, 0;
0, 0, 0, 0, 0.21, 0.42];
tfa4management(P,1,0.9,1.1);% find the minimum absolute value

P_new = [0, 0, 0, 0, 21, 115;
0.32, 0.08, 0, 0, 0, 0;
0, 0.31, 0.11, 0, 0, 0;
0, 0, 0.41, 0.39, 0, 0;
0, 0, 0, 0.16-0.0502, 0.44, 0;
0, 0, 0, 0, 0.21, 0.42];
eigs(P_new,1);

N = [6204; 2335; 823; 482; 119; 37];
%% iteration
Nn = N;
for n = 1:100
    N = P_new*N; 
    Nn = [Nn, N];
end

%% plot
plot(0:100,Nn,'linewidth',2)
xlabel('Time [days]')
ylabel('Population count [individuals]')
legend('Eggs','Neonate larvae','Later-instar larvae','Pupae','Juveniles','Adults')
   