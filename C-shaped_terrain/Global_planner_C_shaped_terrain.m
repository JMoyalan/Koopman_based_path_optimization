%% Data-Driven Optimal navegation for two dimensional single integrator
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontSize',25)
rng(2141444)
%% Parameters Simulation

f = @(t, x, u)([u(1,:);u(2,:)]); %single integrator

A_fin = readmatrix('Terrain Map.xlsx','Sheet','C_shaped');
max_A_fin = max(A_fin,[],"all");
A_fin = A_fin./max_A_fin;
save A_fin A_fin;
tic;
n=2; % dimension
dim = n;

A_fin_rows = size(A_fin,1);
A_fin_columns = size(A_fin,2);
Dom_columns = [-2.5 2.5]; % Domain x1
Dom_rows = [-2.5 2.5]; % Domain x2
dis = 6;    %dis = 12 for vel bounds = 1; dis = 10 for vel bounds = 0.3;
nb_IC =2e4; % number the initial conditions for trining


deltaT = 0.01;
f_d = @(t,x,u) (x + deltaT*f(t, x, u));



f_d1 = @(t,x,u) (x + 1*f(t, x, u));
R_NL =0.01; %globall weight forData=Driven PF Control
n_grid = A_fin_rows;
Nrbf =A_fin_rows*A_fin_columns;   % numer of bases totally

gammma = 1e7; %1e7

%% ************** Phase 1. K approximation  ******************
%  ********* Basis functions **********
sig = (0.49)*(dis/(A_fin_rows - 1)); % sigma in  function of distance
cent_rows = linspace(Dom_rows(1), Dom_rows(2), A_fin_rows);
cent_columns = linspace(Dom_columns(1), Dom_columns(2), A_fin_columns);
[cent_rowsx,cent_columnsy] = ndgrid(cent_rows,cent_columns);
cent = [cent_rowsx(:),cent_columnsy(:)]';



Psi = @(X)(GaussRBF1(X, cent, sig));


Psidot = @(X,dX)(GaussRBF_dot_no_drift(X,dX, cent, sig));
Psi_partial = @(X,pp)(GaussRBF_partial(X, pp, cent, sig));

% Start point
x0 = [-0.875;0.125];
% Goal or Target 
eq_0 = [1.625;0.125];


% [~,idx] = min(vecnorm(cent-eq_0)); % Find the center closest to the equilibrium point
[~,idx] = min(sqrt(sum((cent-eq_0).^2,1))); % Find the center closest to the equilibrium point
cent(:,idx) = eq_0;
[idx] = find(sqrt(sum((cent-eq_0).^2,1))<=(1*sig)); % Find the centers within circle


% index of obstacles equal to one
l = length(idx);
A_Transpose = A_fin';   %To resolve the matrix mismatch issue
b = A_Transpose(:);
b(idx,:) = [];

%%%%%%%%% Hard Obstacle

lg_idx = length(idx);
I = eye(Nrbf-lg_idx);


%X0 = rand(n, nb_IC)*Dom(2)*2 - Dom(2);
X0 = rand(n, nb_IC)*Dom_rows(2)*2 - Dom_rows(2);
[X , Y, Xdot ] = Data(f, X0, deltaT, [0;0]);
[X1 , Y1, X1dot ] = Data(f, X0, deltaT, [1;0]);
[X2 , Y2, X2dot ] = Data(f, X0, deltaT, [0;1]);




PsiX =  Psi(X);
PsiY =  Psi(Y);
PsiX1 =  Psi(X1);
PsiY1 =  Psi(Y1);
PsiX2 =  Psi(X2);
PsiY2 =  Psi(Y2);
data = [X Y X1 Y1 X2 Y2];

[K_f] = K_Operator_positivity(PsiX,PsiY);

[K_fg1] = K_Operator_positivity(PsiX1,PsiY1);

[K_fg2] = K_Operator_positivity(PsiX2,PsiY2);
PF_f = K_f';
PF_fg1 = K_fg1';
PF_fg2 = K_fg2';


d = zeros(Nrbf,1);
qx = @(x1,x2)(x1^2 + x2^2); % cost function associated to the states

for i = 1:Nrbf
    dfunc =@(x1,x2)((x1.^2 + x2.^2).*exp(-(1/sig)^2*((x1-cent(1,i)).^2+(x2-cent(2,i)).^2)));
    d(i) = integral2(dfunc,Dom_rows(1),Dom_rows(2),Dom_columns(1),Dom_columns(2));
end


D = eye(Nrbf);

% ******  Remove rows and colums the origen  ******
D(idx,:) = []; D(:,idx) = [];
PF_f(idx,:) = []; PF_f(:,idx) = [];
PF_fg1(idx,:) = []; PF_fg1(:,idx) = [];
PF_fg2(idx,:) = []; PF_fg2(:,idx) = [];
d(idx) = [];
h = zeros(Nrbf,1);
for iii = 1:size(cent,2)

    p = cent(:,iii);

    if p(1,1) > x0(1,1)-0.5
        if p(1,1) < x0(1,1)+0.5
            if p(2,1) > x0(2,1)-0.5
                if p(2,1) < x0(2,1)+0.5
      h(iii,1) = 1e3; %1e3
                end
            end
        end
    end

end
    
h(idx) = [];
%% ****** Convex  Optimization Problem  ******
%cvx_solver sedumi
tic;
cvx_begin 
variable v(Nrbf-lg_idx,1) %nonnegative
variable w1(Nrbf-lg_idx,1)
variable w2(Nrbf-lg_idx,1)
%cvx_solver sedumi
%variable gammma
ctr_cost = 0;
for i = 1:(Nrbf-lg_idx)
    for j = 1:(Nrbf-lg_idx)
        if i==j
            ctr_cost = ctr_cost + D(i,j)*quad_over_lin(w1(i),v(j)) + D(i,j)*quad_over_lin(w2(i),v(j));
        end
    end
end


minimize(0.001*d'*v + 0.999*b'*v + R_NL*ctr_cost)
% minimize(0.001*d'*v + 0.999*b'*v)

subject to

alpha1 = 1.0;
alpha2 = 0.3;
(I-alpha1*PF_fg1)*w1 + (I-alpha1*PF_fg2)*w2 >= h;
v >= zeros(Nrbf-lg_idx,1);
b'*v <= gammma;
%b1'*v == 0;
abs(w1) <= alpha2*v;
abs(w2) <= alpha2*v;

cvx_end
toc;
cvx_optval

% sum_v = sum(v);
% v = v./sum_v;
% ****  Controller ****
K1t = w1./v;
K2t = w2./v;

K01 = zeros(Nrbf,1);
K01(setdiff(1:end,idx)) = K1t;
K1 = K01;
K02 = zeros(Nrbf,1);
K02(setdiff(1:end,idx)) = K2t;
K2 = K02;
K12 = [K1;K2];

% Rearrangement diemension with the target bases equal to zero
vN = zeros(Nrbf,1); wN1 = zeros(Nrbf,1); wN2 = zeros(Nrbf,1);
vN(setdiff(1:Nrbf, idx)) =v;
wN1(setdiff(1:Nrbf, idx)) =w1;
wN2(setdiff(1:Nrbf, idx)) =w2;


Al = zeros(2,2);
Bl = eye(2);

R_L= 0.1;
Kll =lqr(Al,Bl,eye(n), R_L); %feedback linear gain
Fcl = (Al-Bl*Kll);
P_Lyap = lyap(Fcl,eye(n));
%%
r =50000000;
rho_L=  @(x) max((((x'*P_Lyap*x)^(-3))-r),0);
%K_linear = K_lqr;   % K_local using Bowen 
K_linear = Kll; % K_local using linearized model
%% Simulation
%load exp2_0505
close all

Tmax = 22000;
Nsim = Tmax/1;

rng(2156488)



Ninit = size(x0,2);
% x0_second = 1*rand(n, Ninit) - 0.5;
x0_second = [0;0];
figure(1)
view(225,50)
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$b(\textbf{x})$','interpreter','latex');
set(gca,'fontsize',40)
hold on;
C0 = zeros(size(A_fin,1),size(A_fin,2),3);
for i = 1:size(A_fin,1)
    for j = 1:size(A_fin,2)

        if A_fin(i,j) == 0
            C0(i,j,1) = 0.9290;
            C0(i,j,2) = 0.6940;
            C0(i,j,3) = 0.1250;
        end
        if A_fin(i,j) >= 0.05 && A_fin(i,j) < 0.1
            C0(i,j,1) = 0.3;
            C0(i,j,2) = 0.5;
            C0(i,j,3) = 0.1;
        end
        if A_fin(i,j) >= 0.1 && A_fin(i,j) < 0.2
            C0(i,j,1) = 0.2;
            C0(i,j,2) = 0.3;
            C0(i,j,3) = 0.05;
        end
        if A_fin(i,j) >= 0.2 && A_fin(i,j) < 0.3
            C0(i,j,1) = 0.4;
            C0(i,j,2) = 0.15;
            C0(i,j,3) = 0.045;
        end
        if A_fin(i,j) >= 0.3 && A_fin(i,j) < 0.4
            C0(i,j,1) = 0.2;
            C0(i,j,2) = 0.075;
            C0(i,j,3) = 0.0225;
        end
        if A_fin(i,j) >= 0.4 && A_fin(i,j) < 0.5
            C0(i,j,1) = 0.1;
            C0(i,j,2) = 0.0375;
            C0(i,j,3) = 0.01125;
        end
        if A_fin(i,j) >= 0.5 && A_fin(i,j) < 0.6
            C0(i,j,1) = 0.25;
            C0(i,j,2) = 0;
            C0(i,j,3) = 0;
        end
        if A_fin(i,j) >= 0.6 && A_fin(i,j) < 0.7
            C0(i,j,1) = 0.5;
            C0(i,j,2) = 0;
            C0(i,j,3) = 0;
        end
        if A_fin(i,j) >= 0.7 && A_fin(i,j) < 0.8
            C0(i,j,1) = 0.75;
            C0(i,j,2) = 0;
            C0(i,j,3) = 0;
        end
        if A_fin(i,j) >= 0.8 && A_fin(i,j) < 0.9
            C0(i,j,1) = 1;
            C0(i,j,2) = 0;
            C0(i,j,3) = 0;
        end
        if A_fin(i,j) >= 0.9
            C0(i,j,1) = 1;
            C0(i,j,2) = 0;
            C0(i,j,3) = 0;
        end
    end
end
surf(cent_columns,cent_rows,A_fin,C0);
alpha 0.5


syms X [2, 1] 


V= (X-[eq_0(1);eq_0(2)])'*P_Lyap*(X-[eq_0(1);eq_0(2)]);

g = matlabFunction(V);
fc = fcontour(g,'LineColor',[0 0.4470 0.7410],'LineWidth',5);
fc.LevelList = (1/r)^(1/3); 

V1 = (X-[x0(1);x0(2)])'*P_Lyap*(X-[x0(1);x0(2)]); %%% terrain 2 postion 1

g1 = matlabFunction(V1);
fc1 = fcontour(g1,'LineColor',[1 0 1],'LineWidth',5);
fc1.LevelList = (1/r)^(1/3);


for j = 1:Ninit
    x_true = x0(:,j);
    x_true_second = x0_second(:,j);
    x_open = x0(:,j);
    u_value = []; 
    u_second_value = [];
    
   for i = 0:Nsim-1

          
       
        u_value(1, i+1) = K1'*Psi(x_true(:,end));
        u_value(2, i+1) = K2'*Psi(x_true(:,end));
     

        x_true = [x_true, f_d(0,x_true(:,end), u_value(:, end))];
        x_true_second = [x_true_second, u_value(:, end)];
        u_second_value(1,i+1) = (x_true_second(1,i+2)-x_true_second(1,i+1))/deltaT;
        u_second_value(2,i+1) = (x_true_second(2,i+2)-x_true_second(2,i+1))/deltaT;

        
        
    end
    u_value(:,i+2) =  0;
    u_second_value(:,i+2) = 0; 
   
    
 
    rr1 = plot(x_true(1,:),x_true(2,:),'k','linewidth',5); hold on
    
end

lgr = legend([fc, fc1], 'Target Set', 'Initial Set','interpreter','latex');

hold off;

xlim([Dom_columns(1) Dom_columns(2)]);
ylim([Dom_rows(1) Dom_rows(2)]);



%% Saving csv files for quadsdk
x_bs = [x_true',x_true_second'];
u_bs = u_second_value';


state_filename = 'states_traj.csv';
ctrl_filename = 'ctrl_traj.csv';

writematrix(x_bs, state_filename);
writematrix(u_bs, ctrl_filename);

toc;


