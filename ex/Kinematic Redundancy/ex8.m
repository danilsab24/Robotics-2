clear all
clc

%% Symbolic setup for a 4R planar robot
syms q1 q2 q3 q4 l

p = [ l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3)+l*cos(q1+q2+q3+q4);
      l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3)+l*sin(q1+q2+q3+q4);];
J = jacobian(p,[q1,q2,q3,q4]);

%% Evaluate J at q1=q2=q3=q4=0, l=0.5
J_subs = double(subs(J,[q1,q2,q3,q4,l],[0,0,0,0,0.5]));

v = [0;10];

%% Naive pseudo-inverse solution
J_subs_pinv = pinv(J_subs);
q_dot_unconstrained = J_subs_pinv * v

%% Define velocity bounds
% problem statement says: |q̇₁| ≤ 4, |q̇₂| ≤ 2, |q̇₃| ≤ 1, |q̇₄| ≤ 1
constraint_min = [-4; -2; -1; -1];
constraint_max = [ 4;  2;  1;  1];

%% Call the SNS routine
[new_values, old_values] = SNS(v, {J_subs}, constraint_min, constraint_max, 0);

disp('==================================')
disp('Old (unconstrained) joint velocities:')
disp(old_values)
disp('Feasible joint velocities (SNS) satisfying constraints:')
disp(new_values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNS function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_values, old_values] = SNS(task, J_array, constraint_min, constraint_max, velocity_initial_condition)
    % This algorithm can handle velocity or acceleration (with dJ).
    %
    %  Inputs:
    %    - task: desired velocity or acceleration
    %    - J_array: {J} (for velocity) or {J,dJ} (for acceleration)
    %    - constraint_min: vector of lower bounds for each joint
    %    - constraint_max: vector of upper bounds for each joint
    %    - velocity_initial_condition: (used only if acceleration)
    %
    %  Outputs:
    %    - new_values: feasible joint velocities/accelerations
    %    - old_values: original unconstrained solution from pinv(J)*task
    
    n_joint = length(constraint_min);
    task    = reshape(task, [], 1);  % ensure column vector
    
    constraint_min = reshape(constraint_min, n_joint, 1);
    constraint_max = reshape(constraint_max, n_joint, 1);
    
    init_value  = -999;
    new_values  = init_value*ones(n_joint,1);

    if length(J_array) == 1
        % Velocity mode
        disp("---- SNS for velocity")
        J         = J_array{1};
        type_task = 'v';
        old_values = pinv(J)*task;
    elseif length(J_array) == 2
        % Acceleration mode
        disp("---- SNS for acceleration")
        J         = J_array{1};
        dJ        = J_array{2};
        type_task = 'a';
        velocity_initial_condition = reshape(velocity_initial_condition, n_joint, 1);
        old_values = pinv(J)*(task - dJ*velocity_initial_condition);
    else
        disp("Error: J_array must be {J} or {J, dJ}")
        return
    end

    for i = 1 : n_joint
        fprintf("\n==========> Loop number: %d\n", i)
        
        % Each iteration, compute the unconstrained solution for the
        % *current* J (which has some columns removed in previous loops).
        if type_task == 'v'
            result = pinv(J)*task;
        else
            % type_task == 'a'
            result = pinv(J)*(task - dJ*velocity_initial_condition);
        end
        
        % result might shrink each iteration, so just ensure column shape:
        result = reshape(result, length(result), 1);
        disp("Result: ")
        disp(result)
        
        % Check for any joint that violates constraints
        violating_max_logical = (result > constraint_max);
        violating_min_logical = (result < constraint_min);
        
        % Max violation distance
        violating_max_value = max( (result - constraint_max) .* violating_max_logical );
        % Min violation distance (will be <= 0 if it actually violates)
        violating_min_value = min( (result - constraint_min) .* violating_min_logical );
        
        fprintf("Violating min: %f\n", violating_min_value)
        fprintf("Violating max: %f\n", violating_max_value)
        
        % If no violations, we can finalize
        if all(~violating_min_logical) && all(~violating_max_logical)
            % Fill in the new_values for any joints not yet assigned
            logic_subs = (new_values == init_value);
            % 'result' belongs only to the reduced set of joints (the ones not saturated yet).
            idx_r = find(logic_subs);
            new_values(idx_r) = result; 
            break
        end
        
        % Pick whichever violation is "largest" in absolute value.
        if abs(violating_min_value) > abs(violating_max_value)
            % Identify which joint index triggered that min violation
            idx_candidates = find( (result - constraint_min) == violating_min_value );
            violating_joint_index = idx_candidates(1); % pick the first
            violation_type = 'violating_min';
        else
            % Identify which joint index triggered that max violation
            idx_candidates = find( (result - constraint_max) == violating_max_value );
            violating_joint_index = idx_candidates(1); % pick the first
            violation_type = 'violating_max';
        end
        
        fprintf("--- Joint %d exceeded bounds ==> %f\n", ...
                violating_joint_index, result(violating_joint_index));
        
        % Saturate that joint to its min or max
        if strcmp(violation_type,'violating_min')
            new_values_index = constraint_min(violating_joint_index);
        else
            new_values_index = constraint_max(violating_joint_index);
        end
        fprintf("--- Joint %d now set at %f\n", violating_joint_index, new_values_index)
        
        % Record the saturated value in new_values
        % But we have to figure out which joint in the "global" sense
        % that this corresponds to.  We can do that by scanning for the
        % next -999 in new_values or by storing an index map. For simplicity
        % here, we’ll just do it with a consistent approach:
        global_indxs = find(new_values == init_value); % unsaturated joints so far
        global_joint = global_indxs(violating_joint_index); 
        new_values(global_joint) = new_values_index;
        
        % Subtract that from the task
        selected_J = J(:,violating_joint_index);
        if type_task == 'v'
            % velocity
            task = task - selected_J * new_values_index;
        else
            % acceleration
            task = task - selected_J * new_values_index;
            dJ(:, violating_joint_index) = [];
            velocity_initial_condition(violating_joint_index) = [];
        end
        
        % Now remove that DOF from J and from constraints
        J(:, violating_joint_index) = [];
        constraint_min(violating_joint_index) = [];
        constraint_max(violating_joint_index) = [];
        
    end % for i=1:n_joint
    
    % If for some reason we used all n_joint loops, the new_values might still
    % have some '-999' entries if not all had to saturate. Fill them with
    % the final partial solution from the last iteration:
    leftover_inds = (new_values == init_value);
    if any(leftover_inds)
        % The last 'result' is for the smaller dimension J
        new_values(leftover_inds) = result; 
    end
    
    new_values = reshape(new_values, n_joint,1);  % final shape
end
