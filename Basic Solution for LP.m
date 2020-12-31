%% Defining the problem
% max 3x_1+5x_2
% subject to x_1         + x_3 = 4
%            2x_2             + x_4 = 12
%            3x_1 + 2x_2            + x_5 = 18

%% Clear workspace
clc
clearvars
close ('all')

%% Defining required matrices from from the prblem details
% X = [x1 x2 x3 x4 x5];
sv = 3;             % no of slack variable
c = [3 5 0 0 0];    % objective function coefficient
b = [4 12 18];      % constraints RHS
sv = 3;             % noof slack variables introduced
A = [1 0;
    0 2 ;
    3 2 ];          % constraints coeffcient
augA = [A eye(sv)];

%% Simplex algorithm implementation

% Finding out the initial slack variable to be considered
for i = 1:size(augA,2)
    temp = max(augA(:,i));              % checking for identity matrix in column
    if temp == 1
        count = 0;
        for j = 1:size(augA,1)          % Confirmaing identitiy matrix
            if (augA(j,i) == 0 || augA(j,i) == 1)
                count = count + 1;
            end
        end
    else
        continue
    end
    if count == size(augA,1)
        RI(i) = i;
    end
end

% Obtaining indices for the slack variables
ino = 1;                                % iteration number under consideration
R(:,ino) = RI(RI~=0);                   % first iteration required slack variables to be considered

% Obtaining the required matrices for the clcultions
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
while max(Rcheck(ino,:)) > 0
    if max(Rcheck(ino,:)) > 0
        tempc = find(Rcheck(ino,:) == max(Rcheck(ino,:)));          % Finding out the new variable entering
        theta(:,ino) = zeros((size(augA,1)),1);
        theta(:,ino) = Uc(:,ino)./augA(:,tempc);
        for i = 1:length(theta(:,ino))
            if (theta(i,ino) < 0 || theta(i,ino) == inf)
                theta(i,ino) = NaN;
            end
        end
        tempr = find(theta(:,ino) == min(theta(:,ino)));            % finding out the basic variable leaving
    end
    
    % record for slack variable values for tracking
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

disp('Values of decision variables and slack variables for optimal solution are:')
for i = 1:length(R(:,end))
    formatSpec = 'X_%d  = %d\n';
    fprintf(formatSpec,R(i,end),Uc(i,end))
end
disp('Other variable have the value as 0.')

formatSpec = '\nOptimal solution for objective function with simplex method observed to be %4.2f\n';
fprintf(formatSpec,z(end))