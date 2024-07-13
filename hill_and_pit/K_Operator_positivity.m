function K = K_Operator_positivity(PsiX, PsiY,  gamma  )
% This function compute the approximaiton of Koopman Operator
% PsiX:= Lited input data
% PsiY:= Lifted output data
% gamma:= Factor of regularization
%
%EDMD
M = size(PsiX,2); K = size(PsiX,1);
G = 0;  A = 0;
%
for i = 1:M
    G = G + PsiX(:,i)*PsiX(:,i).';
    A = A + PsiX(:,i)*PsiY(:,i).';
end
G = G/M;
A = A/M;

% if  nargin == 2
%     if rcond(G) > eps
%         K = G\A;
%     else
        K1 = pinv(G)*A;
%     end

% 
% for ii = 1:size(K1,1)
%     for jj = 1:size(K1,2)
%     if K1(ii,jj) < 0
%         K1(ii,jj) = 0;
%     end
%     end
% end

R = zeros(size(K1,1),1);
K = zeros(size(K1));
for i = 1:size(K1,1)
    R(i,1) = sum(K1(i,1:end));
    K(i,:) = K1(i,:).*(1/R(i,1));
end

%K = K';
% 
% elseif nargin == 3
%     yalmip('clear');
%     %*********  Constriant Least-Square Problem, YALMIP. *****
%     Kt = sdpvar(K, K,'full');
%     Objective =  norm(G*Kt-A, 'fro')+ gamma*(norm(Kt, 'fro'));
%     Constraints = [];
%     opt = sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);
%     optimize(Constraints, Objective, opt)
%     %*********
%     K = value(Kt);
% else
%     disp('Error')
% end
end