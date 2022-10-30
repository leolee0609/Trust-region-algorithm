function x = TrustRegion(f, vars,max_Del, ini_Del, x_ini, rou_1, rou_2, gm_1, gm_2, yt, max_iter)
x = x_ini;
del = min(ini_Del, max_Del);
% solve for the hessian matrix and gradient vector of f
B = hessian(f);
b = gradient(f);

% set the stopping criterion to be that #iterations reaches max_iter
iter = 0;
while iter < max_iter
    b_x = double(subs(b, vars, x'));  % current gradient b
    B_x = double(subs(B, vars, x'));  % current hessian B
    
    d = subproblem(B_x, b_x, del);  % solve the subproblem
    disp(iter);
    disp(x);
    f_x_k = double(subs(f, vars, x'));  % f(x^k)
    f_x_k_1 = double(subs(f, vars, (x + d)'));  % f(x^(k+1))
    
    m_k_0 = f_x_k;  % compute m_k(0)
    m_k_d = f_x_k + b_x' * d + (1/2) * d' * B_x * d;  % compute m_k(d)
    
    % compute the measurement
    rou_k = (f_x_k - f_x_k_1) / (m_k_0 - m_k_d);
    
    % update the radius del
    if rou_k < rou_1
        del = gm_1 * del;
    elseif rou_k > rou_2 && norm(d) == del
        del = min(gm_2 * del, max_Del);
    end
    
    % check if we need to update x
    if rou_k > yt
        x = x + d;
    end
    iter = iter + 1;
end


