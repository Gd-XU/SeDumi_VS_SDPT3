% Set up problem parameters
n_values = [3,5,10,20,30,50,100]; % Increasing matrix sizes
m_values = [1,3,5,10,20,30]; % Constraints of dimension m
A_full = randn(max(n_values),max(n_values),max(m_values));
num_problems = length(n_values) * length(m_values);
solver_names = {'sedumi', 'sdpt3'}; % Solvers to compare
tolcontrol = {'sedumi.eps', 'sdpt3.gaptol'}; % Tolerance control method
solver_times = zeros(length(solver_names), num_problems);
tol = 1e-8; % Tolerance for the Lyapunov function

% Solve problems for each combination of n and m
problem_idx = 1; % Problem index
for i = 1:length(n_values)
    n = n_values(i);
    for j = 1:length(m_values)
        m = m_values(j);
        AQ = A_full(1:n,1:n,1:m);
        P = sdpvar(n); % Symmetric nxn matrix
        constraints = [ -eye(n) + P >= 0 ];
        for q = 1:m
            A = AQ(:,:,q);
            constraints = [constraints, -eye(n) - A'*P - P*A >= 0]; 
        end
        for k = 1:length(solver_names)
            solver_name = solver_names{k};
            tolcontrol_name = tolcontrol{k};
            options = sdpsettings('solver', solver_name, 'verbose', 0, tolcontrol_name, tol); 
            objective = trace(eye(n)*P);
            tic;
            sol = optimize(constraints, objective, options);
            solver_times(k, problem_idx) = toc;
            memory_usage(k, problem_idx) = monitor_memory_whos();
        end
        problem_idx = problem_idx + 1;
    end
end

% Print solver times and memory usage for each problem
for i = 1:num_problems
    fprintf('n=%d, m=%d:\n', n_values(ceil(i/length(m_values))), m_values(mod(i-1,length(m_values))+1));
    for j = 1:length(solver_names)
        fprintf('\n%s: %.4f seconds\n', solver_names{j}, solver_times(j,i));
        fprintf('Memory usage: %f MB\n', memory_usage(j, i));
    end
    fprintf('\n');
end

% =========Function=========

function [ memory_in_use ] = monitor_memory_whos( )
% This function uses the 'whos' command and evaluates inside the base
% workspace and sums up the bytes.  The output is displayed in MB.

mem_elements = evalin('base','whos');
if size(mem_elements,1) > 0

    for z = 1:size(mem_elements,1)
        memory_array(z) = mem_elements(z).bytes;
    end

    memory_in_use = sum(memory_array);
    memory_in_use = memory_in_use/1048576;
else
    memory_in_use = 0;
end
end
