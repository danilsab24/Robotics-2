clc; clear; close all;

%% ============== 1) Symbolic setup =====================================
syms q1 q2 q3 real
syms q1_dot q2_dot q3_dot real
syms l real

% End-effector position p(q)
p = [ l*cos(q1) + l*cos(q1+q2) + l*cos(q1+q2+q3);
      l*sin(q1) + l*sin(q1+q2) + l*sin(q1+q2+q3) ];

% Jacobian wrt (q1,q2,q3)
J = jacobian(p, [q1,q2,q3]);
J = simplify(J);

% Time derivative of J: J_dot = sum_i dJ/dqi * q_i_dot
dJdq1 = diff(J,q1);
dJdq2 = diff(J,q2);
dJdq3 = diff(J,q3);

J_dot_expr = dJdq1*q1_dot + dJdq2*q2_dot + dJdq3*q3_dot;

% n(q,q_dot) = J_dot(q,q_dot) * q_dot
n_expr = J_dot_expr*[q1_dot; q2_dot; q3_dot];

%% ============== 2) Numerics at q1=0, q2=0, q3=pi/2, q1dot=0.8, q2dot=0, q3dot=-0.8
q1_val=0; q2_val=0; q3_val=pi/2;
q1dot_val=0.8; q2dot_val=0; q3dot_val=-0.8;
l_val=0.5;

J_num = double( subs(J, [q1,q2,q3,l], ...
                     [q1_val, q2_val, q3_val, l_val]) );
J_dot_num = double( subs(J_dot_expr, ...
  [q1,q2,q3,q1_dot,q2_dot,q3_dot,l], ...
  [q1_val,q2_val,q3_val,q1dot_val,q2dot_val,q3dot_val,l_val]) );
n_num = double( subs(n_expr, ...
  [q1,q2,q3,q1_dot,q2_dot,q3_dot,l], ...
  [q1_val,q2_val,q3_val,q1dot_val,q2dot_val,q3dot_val,l_val]) );

fprintf('J_num =\n'); disp(J_num);
fprintf('J_dot_num =\n'); disp(J_dot_num);
fprintf('n_num =\n'); disp(n_num);

%% ============== 3) Unconstrained min-norm acceleration ================
p_ddot_des = [2; 1];
q_ddot_unc = pinv(J_num)*(p_ddot_des - n_num);

fprintf('Unconstrained q_ddot =\n');
disp(q_ddot_unc);

%% ============== 4) SNS to enforce bounds ==============================
Q_min = [-10; -10; -2];  % from figure (3rd joint must be >= -2)
Q_max = [  7;  10; 10];

% The "task" = p_ddot_des, J_array = {J_num, J_dot_num}, q_dot_init = [0.8;0;-0.8]
[ q_ddot_sns, q_ddot_sns_unc ] = SNS_Acceleration( ...
    p_ddot_des, {J_num,J_dot_num}, Q_min, Q_max, [q1dot_val; q2dot_val; q3dot_val] );

fprintf('\n===== SNS Acceleration Results =====\n');
fprintf('Inside SNS, unconstrained solution was:\n');
disp(q_ddot_sns_unc);
fprintf('Final SNS-saturated solution:\n');
disp(q_ddot_sns);


%% ============== The SNS function ======================================
function [new_values, old_values] = SNS_Acceleration(task, J_array, q_min, q_max, q_dot_init)
    % SNS_Acceleration  Saturates per-joint accelerations to respect bounds
    %
    % We solve:
    %    q_ddot = pinv(J)*(task - dJ*q_dot_init),
    % subject to q_min <= q_ddot <= q_max.
    %
    % Inputs:
    %    - task:  desired Cartesian accel (2×1)
    %    - J_array:  {J, dJ}, numeric
    %    - q_min, q_max:  each 3×1 for the example
    %    - q_dot_init:    3×1 current velocities (used in n(q, q_dot))
    % Outputs:
    %    - new_values: final feasible q_ddot
    %    - old_values: unconstrained (pinv) solution

    if length(J_array)~=2
        error('SNS_Acceleration expects {J, dJ} in J_array!');
    end
    J  = J_array{1};
    dJ = J_array{2};

    n_joints = length(q_min);

    % 1) "old_values" = unconstrained min-norm solution
    old_values = pinv(J)*(task - dJ*q_dot_init);

    % 2) We start with all joints unsaturated (marked as -999)
    new_values = -999*ones(n_joints,1);

    % 3) Local copies that shrink each iteration
    JJ = J;     % partial J
    dJJ= dJ;    % partial dJ
    t  = task;  % partial task
    qm = q_min; % local min bounds
    qM = q_max; % local max bounds
    qv = q_dot_init;

    for count=1:n_joints
        partial_sol = pinv(JJ)*(t - dJJ*qv);
        partial_sol = partial_sol(:);

        too_big   = (partial_sol >  qM(1:length(partial_sol)));
        too_small = (partial_sol <  qm(1:length(partial_sol)));

        max_vio = max( partial_sol - qM(1:length(partial_sol)) );
        min_vio = min( partial_sol - qm(1:length(partial_sol)) );

        if all(~too_big) && all(~too_small)
            % no violation => we can store partial_sol in new_values
            free_joints = find(new_values==-999);
            new_values(free_joints) = partial_sol;
            break;
        end

        % pick bigger absolute violation
        if abs(min_vio)>abs(max_vio)
            % saturate min
            idx_candidates = find( (partial_sol - qm(1:length(partial_sol))) == min_vio );
            side = 'min';
        else
            % saturate max
            idx_candidates = find( (partial_sol - qM(1:length(partial_sol))) == max_vio );
            side = 'max';
        end
        sat_index_local = idx_candidates(1);

        free_joints = find(new_values==-999);  % which are still unsaturated
        global_idx  = free_joints(sat_index_local);

        if strcmp(side,'min')
            satur_val = qm(sat_index_local);
        else
            satur_val = qM(sat_index_local);
        end
        new_values(global_idx) = satur_val;

        % subtract from the partial task
        selJ = JJ(:, sat_index_local);
        t    = t - selJ*satur_val;

        % remove that DOF from partial system
        JJ(:, sat_index_local)  = [];
        dJJ(:, sat_index_local) = [];
        qv(sat_index_local)     = [];
        qm(sat_index_local)     = [];
        qM(sat_index_local)     = [];
    end

    % If we still have leftover -999, fill from partial_sol
    leftover = find(new_values==-999);
    if ~isempty(leftover)
        for k=1:length(leftover)
            new_values(leftover(k))= partial_sol(k);
        end
    end
end
