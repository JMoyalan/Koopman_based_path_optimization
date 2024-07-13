function Psi = GaussRBF( X,C, sig)
% X: (dim, N)
% C: (dim, Nrbf)

rows = zeros(size(X, 2), 20);
cols = zeros(size(X, 2), 20);
values = zeros(size(X, 2), 20);
for i = 1:size(X, 2)
    [mink_sum_abs, indexes] = mink(sum(abs(X(:, i)-C), 1), 20);
    rows(i, :) = indexes';
    cols(i, :) = i;
    values(i, :) = exp(-sum( (X(:, i) - C(:, indexes)).^2 ,1)/sig^2)';
end
Psi = sparse(rows, cols, values, size(C, 2), size(X, 2)); 

% tic
% Ctemp = C;
% for i = 1:size(C, 2) % for each basis
%     C = repmat( Ctemp(:,i), 1, size(X,2) );
%     psi = exp(-sum( (X - C).^2 ,1)/sig^2);
%     
%     Psi(i,:) = psi;
% end
% toc
