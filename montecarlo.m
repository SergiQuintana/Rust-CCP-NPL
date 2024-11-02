clear all

close all
clc


% Set Parameters
rng(1)
beta=0.999;
thetas = [9.74, 2.69];
trans = [0.29,0.7,0.01];

states = linspace(1,90,90)';

 montecarlo_simulation(beta,thetas,trans,states)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .

function montecarlo_simulation(beta,thetas,trans,states)

% simulate data for different sample sizes

disp('Solving for n=100')
n = 100;
t = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,t);
data = data(1:end-1,:);
% Estimate full solutions
tic;
[theta100,std100] =get_estimates(beta,choices,data,states,trans);
time100 = toc;
% Estimate ccps
tic;
[thetaccp100,stdccp100] = get_estimates_aguirre(beta,choices,data,states,trans,1);
timeccp100 = toc;

disp('Solving for n=1000')
n = 1000;
t = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,t);
data = data(1:end-1,:);
% Estimate full solutions
tic;
[theta1000,std1000] =get_estimates(beta,choices,data,states,trans);
time1000 = toc;
% Estimate ccps
tic;
[thetaccp1000,stdccp1000] = get_estimates_aguirre(beta,choices,data,states,trans,1);
timeccp1000 = toc;


disp('Solving for n=5000')
n = 5000;
t = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,t);
data = data(1:end-1,:);
% Estimate full solutions
tic;
[theta5000,std5000] =get_estimates(beta,choices,data,states,trans);
time5000 = toc;
% Estimate ccps
tic;
[thetaccp5000,stdccp5000] = get_estimates_aguirre(beta,choices,data,states,trans,1);
timeccp5000 = toc;

disp('Solving for n=10000')
n = 10000;
t = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,t);
data = data(1:end-1,:);
% Estimate full solutions
tic;
[theta10000,std10000] =get_estimates(beta,choices,data,states,trans);
time10000 = toc;
% Estimate ccps
tic;
[thetaccp10000,stdccp10000] = get_estimates_aguirre(beta,choices,data,states,trans,1);
timeccp10000 = toc;

disp('Solving for n=25000')
n = 25000;
t = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,t);
data = data(1:end-1,:);
% Estimate full solutions
tic;
[theta25000,std25000] =get_estimates(beta,choices,data,states,trans);
time25000 = toc;
% Estimate ccps
tic;
[thetaccp25000,stdccp25000] = get_estimates_aguirre(beta,choices,data,states,trans,1);
timeccp25000 = toc;


% Generate the table

output1 = zeros(10,1);
output2 = zeros(10,1);

output1(1,1) = theta100(1);
output1(2,1) = thetaccp100(1);
output1(3,1) = theta1000(1);
output1(4,1) = thetaccp1000(1);
output1(5,1) = theta5000(1);
output1(6,1) = thetaccp5000(1);
output1(7,1) = theta10000(1);
output1(8,1) = thetaccp10000(1);
output1(9,1) = theta25000(1);
output1(10,1) = thetaccp25000(1);

output1(1,2) = std100(1);
output1(2,2) = stdccp100(1);
output1(3,2) = std1000(1);
output1(4,2) = stdccp1000(1);
output1(5,2) = std5000(1);
output1(6,2) = stdccp5000(1);
output1(7,2) = std10000(1);
output1(8,2) = stdccp10000(1);
output1(9,2) = std25000(1);
output1(10,2) = stdccp25000(1);


output2(1,1) = theta100(2);
output2(2,1) = thetaccp100(2);
output2(3,1) = theta1000(2);
output2(4,1) = thetaccp1000(2);
output2(5,1) = theta5000(2);
output2(6,1) = thetaccp5000(2);
output2(7,1) = theta10000(2);
output2(8,1) = thetaccp10000(2);
output2(9,1) = theta25000(2);
output2(10,1) = thetaccp25000(2);

output2(1,2) = std100(2);
output2(2,2) = stdccp100(2);
output2(3,2) = std1000(2);
output2(4,2) = stdccp1000(2);
output2(5,2) = std5000(2);
output2(6,2) = stdccp5000(2);
output2(7,2) = std10000(2);
output2(8,2) = stdccp10000(2);
output2(9,2) = std25000(2);
output2(10,2) = stdccp25000(2);

output1 = round(output1,3);
output2 = round(output2,3);

time = zeros(10,1);
time(1) = time100;
time(2) = timeccp100;
time(3) = time1000;
time(4) = timeccp1000;
time(5) = time5000;
time(6) = timeccp5000;
time(7) = time10000;
time(8) = timeccp10000;
time(9) = time25000;
time(10) = timeccp25000;

observations = zeros(10,1);
observations(1) = 100;
observations(2) = 100;
observations(3) = 1000;
observations(4) = 1000;
observations(5) = 5000;
observations(6) = 5000;
observations(7) = 10000;
observations(8) = 10000;
observations(9) = 25000;
observations(10) = 25000;

estimation = zeros(10,1);
estimation(2) = 1;
estimation(4) = 1;
estimation(6) = 1;
estimation(8) = 1;
estimation(10) = 1;

names = ["Observations";"Estimation";"Theta R";"stdR";"Theta 1";"std1";"Time"];

time = round(time,2);

output = table(observations,estimation,output1(:,1),output1(:,2),output2(:,1),output2(:,2),time,'VariableNames',names)
saveTableToLatex(output, 'Output/montecarlo_table.tex')


end



function [thetahat,stdtheta] = get_estimates_aguirre(beta,choices,data,states,trans,iterations)


% Get transition matrix
trans_matrix = get_trans_matrix(trans);

% estimate the model
[thetahat,hessianhat] =aguirremira(data,choices,states,trans_matrix,beta,iterations);

% get st errors
stdtheta=sqrt(diag(inv(hessianhat)));

end

function [thetahat,stdtheta] = get_estimates(beta,choices,data,states,trans)


% compute transition frequencies
choices(choices==0) =2;
trans_matrix = get_trans_matrix(trans);
theta1 =[0,0];

% estimate the model

[thetahat,fval13,exti,output,grad,hessianhat]= fminunc(@(thetas) loglike_full(thetas,beta,trans_matrix,choices,data,states), theta1);

% get st errors
stdtheta=sqrt(diag(inv(hessianhat)));

end

function out = loglike_full(thetas,beta,trans_matrix,choices,data,states)

% Find the fixed point
EV = valuefunction(states,trans_matrix,beta,thetas);

[T,N] = size(data);

% loop over all periods
ft =zeros(T,1);

for t=1:T
    
    % get the conditional value functions
    
    vjt = get_conditional(data(t,:),thetas,EV,beta,trans_matrix);
    
    % set non-replacement as the base category: 
    
    vjt  = vjt - vjt(:,2);
    
    % get the ccps    
    ccps = exp(vjt)./(exp(vjt(:,1))+exp(vjt(:,2)));
    
    % use the ccp of the corresponding choice
    loglike  =  zeros(N,1);
    for i=1:N
    
        loglike(i) = log(ccps(i,choices(t,i)));
        
    end
    ft(t) = sum(loglike);
end
    out = -sum(ft);
end



function vjt = get_conditional(data,theta,EV,beta,trans_matrix)
%This function draws a gumbell and perfomrs the chocie
data = data';
% First get utilities
u = get_utilities(theta,data,EV,beta);

C = get_continuation(data,EV,trans_matrix);

% get conditional value function
v1 = u(:,1) + beta*EV(1);
v2 = u(:,2) + beta*C;
vjt = [v1,v2];

end



function trans_matrix = get_trans_matrix(trans)
% This function copmutes the transition probability matrix given the
% transition probabilities

trans_matrix = zeros(90,90);

for i=1:88
    trans_matrix(i,i:i+2) = trans;
end

trans_matrix(89,89) = trans(1);
trans_matrix(89,90) = 1- trans(1);
trans_matrix(90,90) = 1;

end


function EV=valuefunction(states,F,beta,thetas)
% valuefunction.m iterates on the Expected value function to find a fixed point

% Inputs:
% - state - vector of K possible state, e.g [1,90] or [1,175]
% - F - (K x K) transition matrix  
% - beta - discount factor (two possibile calibrations: 0.9999 or 0)
% - theta1:  is a (2 x 1) vector of parameters we want to estimate:
%         - theta(1) is the the replecement cost
%         - theta(2) is the ordinary mantainance cost parameter

% Outputs:
% - out - (K x 1) vector of expected value functions for each of the K states

% In this example, we consider a linear cost function i.e.:
% c= 0.001*theta(2)*x

% size of grid point
k=size(states,1);

% Compute the mantainance costs for each state;
c= 0.001*thetas(2)*(states);

% Tolerance level in value function iteration
toler=1e-10;

% Convergence norm
rule=1;

% Starting guess
EV0=zeros(k,1);

%Iteration
while rule >toler 
expected = F*EV0;
exp1= exp(-c + beta*expected);

exp2= exp(-thetas(1) -c(1) + beta*EV0(1) );

EV1  = log(exp1 + exp2);
   
rule=norm(EV1-EV0) ;

EV0 = EV1;

end

EV=EV1;
end


  
function u =get_utilities(thetas,data,EV,beta)
% This function gets the utilities of each alternative given the
% paramenters and the value functions. 
u_replace = -thetas(1)*ones(length(data),1);
u_notreplace = -0.001*thetas(2)*(data);
u = [u_replace,u_notreplace];
end


function C = get_continuation(data,EV,trans_matrix)

expected = trans_matrix*EV;
C  = ones(length(data),1);


for i=1:length(data)
C(i) = expected(data(i));
end

end


function choice = make_choice(data,theta,EV,beta,trans_matrix)
%This function draws a gumbell and perfomrs the chocie

% First get utilities
u = get_utilities(theta,data,EV,beta);

%Now draw type-1 extreme value shock
e = evrnd(0,1,[length(data),2]);

C = get_continuation(data,EV,trans_matrix);

% get conditional value function
v1 = u(:,1) + beta*EV(1) + e(:,1);
v2 = u(:,2) + beta*C + e(:,2);
v = [v1,v2];

%find the maximum as the choice

[val,choice] = max(v');
choice(choice==2) = 0;
end

function newstate = move_stocastically(data,trans)
if data <89
    r = mnrnd(1,trans,1);
    increase = find(r==1) -1;
    newstate = data+increase;
end
if data == 89
    p = [trans(1),1-trans(1)];
    r = mnrnd(1,p,1);
    increase = find(r==1) -1;
    newstate = data+increase;
end

if data == 90
    newstate = 90;
end
end

function datanew = move_state(choice,data,trans)
%this function moves the state space given the choice:
datanew = zeros(length(data),1);

for i=1:length(data)
    if choice(i) == 1
        datanew(i) = 1;
    else
        datanew(i) = move_stocastically(data(i),trans);
    end

end

end



function [data,choices] = simulate_data(states,theta,trans,beta,n,T)
%this function simulates n buses for t 
tr = get_trans_matrix(trans);
EV=valuefunction(states,tr,beta,theta);
data = ones(T,n);
choices = zeros(T,n);
data(1,:) = randi([1 90],1,n);
for t=1:T
    choices(t,:) = make_choice(data(t,:)',theta,EV,beta,tr);
    data(t+1,:) = move_state(choices(t,:),data(t,:),trans);
end
save('rust_data',"data","choices")
end

function [thetahat,hessianhat] = aguirremira(data,choices,states,trans_matrix,beta,iterations)
% this function takes as input the empirical ccps and performs the
% aguirregabira mira estimation for some predetermined amount of
% iterations.

theta0 = [0,0];

for it=1:iterations

    % Get the current ccps

    if it == 1
        newccp = get_ccps(data,choices);
    end
    if it>1
        newccp = update_ccps(theta0,states,beta,trans_matrix,oldccp);
    end

    % find the optimal values using those ccps

    [thetahat,fval13,exti,output,grad,hessianhat]= fminunc(@(thetas) loglike_ccp(thetas,beta,choices,data,newccp,trans_matrix), theta0);

    theta0 = thetahat
    oldccp = newccp;


end


end

% let's construct the likelihood

function out = loglike_ccp(thetas,beta,choices,data,ccps,trans_matrix)


[T,N] = size(data);

% loop over all periods
ft =zeros(T,1);

% get expected ccps:

expected = trans_matrix*log(ccps);

for t=1:T
    
    loglike  =  zeros(N,1);
    for i=1:N
        % get the payoff difference using ccps so that continuation value
        % cancels out. 
        diff  = -thetas(1) +beta*(-thetas(1) -log(ccps(1))) - (-0.001*thetas(2)*(data(t,i))+beta*(-thetas(1) -expected(data(t,i))));
        ccp_norepalce = 1/(1+exp(diff));
        % construct the likelihood
        loglike(i) =  choices(t,i)*log(1-ccp_norepalce)...
                      +(1-choices(t,i))*log(ccp_norepalce);
        
    end
    ft(t) = sum(loglike);
end
    out = -sum(ft);
end



function ccp = get_ccps(data,choices)
% This function computes the empirical frequencies of the data
total = tabulate(reshape(data.',1,[]));
total = total(:,1:2);
datareplaced = data.*choices;
total_replaced =  tabulate(reshape(datareplaced.',1,[]));
total_replaced = total_replaced(2:end,1:2);
ccp = zeros(90,1);
for i=1:90
    num = total_replaced(total_replaced(:,1)==i,2);
    if isempty(num)
        num = 0;
    end
    den = total(total(:,1)==i,2);
    if isempty(den)
        den = 1;
    end
    ccp(i) = num / den;
    if ccp(i) == 0
        ccp(i) = 0.0001; % replace by a very small number
    end
end

end

function newccp = update_ccps(theta,states,beta,trans_matrix,oldccp)

% get expected ccps

expected = trans_matrix*log(oldccp);

% get payoff difference

diff = -theta(1) + beta*(-theta(1)-log(oldccp(1))) - (-0.001*theta(2)*states + beta*(-theta(1) - expected));

newccp = 1-  (1./(1+exp(diff)));

end


function realccp = get_realccps(thetareal,states,beta,trans_matrix)

% get continuation value
EV=valuefunction(states,trans_matrix,beta,thetareal);

% First get utilities
u = get_utilities(thetareal,states);

C = get_continuation(states,EV,trans_matrix);

% get conditional value function
v1 = u(:,1) + beta*EV(1);
v2 = u(:,2) + beta*C;

realccp = 1-  (1./(1+exp(v1-v2)));

end

function saveTableToLatex(T, filename)
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write the beginning of the tabular environment
    fprintf(fid, '\\begin{tabular}{%s}\n', repmat('c', 1, width(T)));
    fprintf(fid, '\\hline\n');
    
    % Write the table column names
    colNames = T.Properties.VariableNames;
    fprintf(fid, '%s \\\\\n', strjoin(colNames, ' & '));
    fprintf(fid, '\\hline\n');
    
    % Write the table data
    for i = 1:height(T)
        row = table2cell(T(i, :)); % Convert the row to a cell array
        rowStr = cell(size(row));
        for j = 1:length(row)
            if isnumeric(row{j})
                rowStr{j} = num2str(row{j});
            else
                rowStr{j} = row{j};
            end
        end
        fprintf(fid, '%s \\\\\n', strjoin(rowStr, ' & '));
    end
    
    % Write the end of the tabular environment
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    
    % Close the file
    fclose(fid);
end


