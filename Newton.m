function x = Newton(f, vars, max_iter, x_ini)
x = x_ini;

H = hessian(f);  % get the Hessian matrix of f
g = gradient(f);  % get the gradient vector of f

% start iteration until reaches max_iter
iter = 0;
while iter < max_iter
    disp(iter);
    disp(x);
    % get the current Hessian matrix and gradient vector of f
    H_x = double(subs(H, vars, x'));
    g_x = double(subs(g, vars, x'));

    % get the inverse of H to get d
    d = -inv(H_x) * g_x;

    % update x
    x = x + d;

    iter = iter + 1;
end

end
