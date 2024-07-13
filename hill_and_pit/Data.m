function [X, Y, Xdot] = Data(VecField, X0, Delta_t, u, solver)
% This fucniton use get the samples of one-step
%Euler method ofr obtain data
if nargin == 4
    solver =  'Euler';
end
switch solver
    case 'Euler'
        f = @(x) (x + Delta_t*VecField(0, x, u));
        fdot = @(x) (VecField(0, x, u));
        X = X0;
        Y = f(X) ;
        Xdot = fdot(X);
        
    case  'Ode45'
         X = [];  Y = [];
           tspan = [0  Delta_t];
           options = odeset('RelTol',1e-9,'AbsTol',1e-300);
        for i = 1:size(X0, 2)          
            [t, z]  = ode45(@(t, x)VecField(t, x, u), tspan, X0(:, i), options);
            X = [X z(1:end-1,:)'];
            Y = [Y z(2:end,:)'];
        end
end
end