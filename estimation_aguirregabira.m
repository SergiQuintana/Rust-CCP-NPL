
% Code to simulate rust data

clear all
close all
clc


% define some parameters

beta=0.999;
trans = [0.29,0.7,0.01];
states = linspace(1,90,90)';

% load the data

load('rust_data.mat');
[T,N] = size(data);
data = data(1:T-1,:);
trans_matrix = get_trans_matrix(trans);
thetareal = [12, 2];

iterations = 1000;
aguirremira(data,choices,states,trans_matrix,beta,iterations)


function thetahat = aguirremira(data,choices,states,trans_matrix,beta,iterations)
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

    [thetahat]= fminunc(@(thetas) loglike_ccp(thetas,beta,choices,data,newccp,trans_matrix), theta0);

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

  
function u =get_utilities(thetas,data)
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




