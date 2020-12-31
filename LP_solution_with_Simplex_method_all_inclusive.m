%% Defining the problem
% max x_1 + x_2 + x_3 + x_4
% subject to x_1 + x_2      + x_5 = 3
%            x_3 + x_4      + x_6 = 2

%% Clear workspace
clc
clearvars
close ('all')

%% Defining required matrices from from the prblem details
% X = [x1 x2 x3 x4 x5 x6];
sv = 2;                 % no of slack variable
c = [1 1 1 1 0 0];      % objective function coefficient
b = [3 2];              % constraints RHS

A = [1 1 0 0;
    0 0 1 1] ;          % constraints coeffcient
augA = [A eye(sv)];
set = 1;                % termination criteria for feasible solution 
ino_opt = 0;            % iterations for alternate optimum 
breakout = 0;           % parameter for termination of the alternate optimum loop

%% Simplex algorithm implementation
% Finding out the initial slack variable to be considered
ino = 1;                                % iteration number under consideration
R = [(size(A,2)+1):(size(augA,2))]';    % initial slack variables to be considered

% Obtaining the required matrices for the calculations
for i = 1:length(R)
    augC(i) = c(R(i));
    augB(:,i) = augA(:,R(i));
end

% Calculation of zj for the first iteration
z(ino,:) = augC*inv(augB)*augA;           % Calculation of zj for first iteraion
Uc(:,ino) = inv(augB)*b';                 % Calculaion of updated constraints at first iteration

z(ino,end+1) = augC*inv(augB)*b';         % Calculation of the cost function at the interval
Rcheck(ino,:) = c - z(ino,1:end-1);       % Defining the matix for the ietraion decision for optimal soluation
recaugA = [ones(1,size(augA,2))*ino; augA]; % matrix for record of change of coeff in consraints for iteractions
recaugB = [ones(1,size(augB,2))*ino; augB]; % matrix for record of the chnage in values of RHS for the constraints in each iterations

%% Defining second onwards iterations for the simplex tabular form method
while max(Rcheck(ino,:)) >= 0
    % Checking for the infeasible solution and alternate optimum
    if max(Rcheck(end,:)) == 0
        base_var = (1:size(A,2));
        slack_var = (size(A,2)+1:size(augA,2));
        temp_ix = (1+size(augA,2))*ones(1,size(base_var,2));
        % Confirming for slack variable requirement in final optimal
        % solution possible only if (no of sv > no of base variable)
        if sv <= size(base_var,2)
             for i = 1:length(base_var)
                if isempty([find(base_var(i) == R(:,end))]) ~= 1 
                    temp_ix(i) = base_var(i);
                elseif ((Rcheck(end,base_var(i)) ~= 0) && (Rcheck(end,base_var(i)) < 0))
                    set = 0;            % Marking condition for infeasible or nonfeasible solution
                    break
                else
                    temp_ix(i) = 0;     % Marking for alternate optimum solution
                end
            end
        end
        % Decison for non feasible and alternate optimum calculations
        if isempty([find(temp_ix == 0)]) == 1
            break
        else
            ino_opt = ino_opt + 1;              % Starting iteraions for alternate optimum
            ne_opt = 1;                         % to move fornext set of optimal solution
            R_opt(:,ino_opt) = R(:,ino);        % Optimum solution variable indexes
            z_opt(:,ino_opt) = z(ino,end);      % Optimum solution valuefor check
            Uc_opt(:,ino_opt) = Uc(:,ino);      % Optimum solution vriable values
            % Using rule of in case of zero entering in simplex after one
            % iteration solution repeats itself hence using it as the
            % criteria formoving forward and updating the solution
            if size(R_opt,2)>2
                for i = 1:(size(R_opt,2)-1)
                    t_opt = R_opt(:,ino_opt) - R_opt(:,i);  % Checking for repeated solution
                    if size((find(t_opt == 0)),1) == sv
                        ino = ino - 1;                      % Removing last iteration producing earlier solution from main loop
                        ino_opt = ino_opt - 1;              % Removing last iteration producing earlier solution from alternate optimum loop
                        ne_opt = 2;                         % introducing variable for entering base variable which has notbeen entered
                        breakout = breakout + 1;            % Variable to avoid continuous loop in case of alternate optimum 
                            if breakout > 2
                                ne_opt = 3;                 % variable to break out the alternate optimum continuous loop 
                                break
                            end
                        break                               % to come out of the for loop at first same value of optimum solution
                    end
                end
            end
            if ne_opt == 2
                for i = 1:length(base_var)
                    if isempty(find(R_opt(:) == i))
                        if Rcheck(end,i) >= 0               % verifying base varible notentered have either zero or postive value while entering       
                            tempc = i;                      % finding base variable not presentin earlier 
                            break
                        else
                            break
                        end
                    end
                end
            elseif ne_opt == 3
                break                                       % termination of the continuous alternative optimum loop
            else
               tempc = find((temp_ix == 0),1);              % new value for entering in case of alternate optimum in normal condition
            end
            % theta calculation to determine the leaving varibable
            theta(:,ino) = zeros((size(augA,1)),1);
            theta(:,ino) = Uc(:,ino)./augA(:,tempc);
            for i = 1:length(theta(:,ino))
                if (theta(i,ino) < 0 || theta(i,ino) == inf)
                    theta(i,ino) = NaN;
                end
            end
            tempr = find(theta(:,ino) == min(theta(:,ino)));
            set = 2;
        end       
    end
    % Calculations for feasible and unbounded solutions
    if max(Rcheck(ino,:)) > 0
        if max(Rcheck(ino,:)) > 0
            tempc = find((Rcheck(ino,:) == max(Rcheck(ino,:))),1);          % Finding out the new variable entering
            theta(:,ino) = zeros((size(augA,1)),1);
            theta(:,ino) = Uc(:,ino)./augA(:,tempc);
            for i = 1:length(theta(:,ino))
                if (theta(i,ino) < 0 || theta(i,ino) == inf)
                    theta(i,ino) = NaN;
                end
            end
            tempr = find(theta(:,ino) == min(theta(:,ino)));            % finding out the basic variable leaving
        end

        % providing condition for the theta check for unboundedness
        if ((isnan(max(theta(:,ino)))) == 1) && ((isnan(min(theta(:,ino)))) == 1)
            set = inf;                  % unbounded solution
            disp('Function to maximize is unbounded ')
            break
        end
        set = 1;
    end
        ino = ino + 1;
        recaugA = [recaugA; ino*ones(1,size(augA,2)); augA];
        recaugB = [recaugB; ino*ones(1,size(augB,2)); augB];

        % Updating basic variables for the solution
        for i = 1:size(R,1)
            if i == tempr
                R(i,ino) = tempc;
            else
                R(i,ino) = R(i,ino-1);
            end 
        end

        Uc(tempr,ino) = Uc(tempr,ino-1)./augA(tempr,tempc);
        augA(tempr,:) = augA(tempr,:)./augA(tempr,tempc);
        for i = 1:size(augA,1)
            if (i ~= tempr)
                Uc(i,ino) = Uc(i,ino-1) - (Uc(tempr,ino)*augA(i,tempc));
                augA(i,:) = augA(i,:) - (augA(tempr,:)*augA(i,tempc));
            end
        end

        for i = 1:size(R,1)
            augC(ino,i) = c(R(i,ino));
            augB(:,i) = augA(:,R(i,ino));
        end

        % Calculation of zj for iteration under consideration
        z(ino,1:end-1) = augC(ino,:)*inv(augB)*augA;            % Calculation of zj for first iteraion
        z(ino,end) = augC(ino,:)*inv(augB)*Uc(:,ino);           % Calculation of the cost function for the iteration
        Rcheck(ino,:) = c - z(ino,1:end-1);
end

if set == inf
    formatSpec = '\nGiven objective function does not have optimal solution for maximization\n';
    fprintf(formatSpec,z(end))
elseif set == 0
    formatSpec = '\nGiven constriants for objective function have infeasible region\n';
    fprintf(formatSpec,z(end))
    formatSpec = '\nGiven objective function have infeasible solution for maximization\n';
    fprintf(formatSpec,z(end))
elseif set == 2
    if z_opt(:) == z_opt(end)
        formatSpec = 'Solution have mutliple optimal solution with optimal function value as %4.2f\n';
        fprintf(formatSpec,z_opt(end))
    end
    formatSpec = '\nMultiple optimal solutions are as below-\n';
    fprintf(formatSpec)
    disp('Values of decision variables and slack variables for optimal solution are:')
    for i = 1:(size(R_opt,2)-1)
        formatSpec = 'Solution (%d) is\n';
        fprintf(formatSpec,i)
        for j = 1:length(R_opt(:,i))
            formatSpec = 'X_%d  = %d\n';
            fprintf(formatSpec,R_opt(j,i),Uc_opt(j,i))
        end
        disp('Other variable have the value as 0.')
    end
else
    disp('Values of decision variables and slack variables for optimal solution are:')
    for i = 1:length(R(:,end))
        formatSpec = 'X_%d  = %d\n';
        fprintf(formatSpec,R(i,end),Uc(i,end))
    end
    disp('Other variable have the value as 0.')
    formatSpec = '\nOptimal solution for objective function with simplex method observed to be %4.2f\n';
    fprintf(formatSpec,z(end))
end