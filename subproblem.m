function d = subproblem(B, b, del)
% first, do eigendecomposition of B
lambd = eig(B);  % get all eigenvalues
lambd1 = min(lambd);  % get the minimum eigenvalue
% see where B is pd by checking its eigenvalues
if lambd1 > 0
    % B is pd and d has norm less than delta, then B has inverse
    d = -inv(B) * b;  % calculate d
    if norm(b) < del
        % d has norm less than delta, then d is optimal, and we are done
        return;
    end
end
% get the eigenvector v1 of lambda1
[V,D] = eig(B);
n = size(B);
n = n(1);  % n is the size of B
% initialize a matrix to store ALL eigenvectors that has eigenvalue
%   lambda1
Lbd1eigv = zeros(n, n);
for i = 1: n
     if D(i, i) == lambd1
         Lbd1eigv(:, i) = V(:, i); 
     end
end
% see if ALL such eigenvectors have dot product 0 with b
check_vec = (b' * Lbd1eigv)';
if norm(check_vec) == 0
    disp(-lambd1);
    % ALL such eigenvectors have dot product 0 with b
    % calculate -d(-lambda1)
    d_lbd1 = 0.0;
    for j = 1: n
        if D(j, j) ~= lambd1
             d_lbd1 = d_lbd1 + ((V(:, j)' * b) / (D(j, j) - lambd1)) * V(:, i);
        end
    end
    disp(norm(d_lbd1))
    if norm(d_lbd1) <= del
        
        tau = sqrt(del.^2 - norm(d_lbd1).^2);
        d = -d_lbd1 + tau * V(:, 1);  % d is the optimal solution
        return;
    end
end
% calculate lower and upper bound
lower_sig_term = 0.0;
for j = 1: n
    if D(j, j) == lambd1
         
         lower_sig_term = lower_sig_term + (V(:, j)' * b).^2;
    end
end
lower_sig_term = lower_sig_term.^(1/2);
lbd_lower = -lambd1 + (1 / del) * lower_sig_term;

lbd_upper = norm(b) / del - lambd1;

lbd_lower = max(lbd_lower, 0);

% implement newton's method to find root
err = 0.0000001;  % the max error of the root
lbd_st = (lbd_upper + lbd_lower) / 2;  % initial guess is between upper and lower bound
lbd_st = lbd_upper;
Idtt = eye(n);
dist_to_x = norm(inv(B + lbd_st * Idtt) * b) - del;  % distance to the x-axis
while abs(dist_to_x) > err
    % get the derivative of the function
    deri = - (b' * inv(B + lbd_st * Idtt)' * inv(B + lbd_st * Idtt).^2 * b) / norm(inv(B + lbd_st * Idtt) * b);
    % update lambda^*
    lbd_st = lbd_st - dist_to_x / deri;
    
    % update distance to the x-axibs
    dist_to_x = norm(inv(B + lbd_st * Idtt) * b) - del;
end
disp(lbd_st);
% get the optimal d
d = -inv(B + lbd_st * Idtt) * b;
return
end