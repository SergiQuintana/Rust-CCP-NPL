
% Code to simulate rust data

clear all
close all
clc


%% 
% Set Parameters

beta = 0.99; 
thetas = [9.74, 2.69];
trans = [0.29,0.7,0.01];

states = linspace(0,445000,90)';
data = linspace(1,50,50);

tr = get_trans_matrix(trans);

EV=valuefunction(states,tr,beta,thetas);
u =get_utilities(thetas,data,states,EV,beta);

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
exp1= exp(-c + beta*EV0);

exp2= exp(-thetas(1) -c(1) + beta*EV0(1) );

EV1  = F * log(exp1 + exp2);
   
rule=norm(EV1-EV0) ;

EV0 = EV1;

end

EV=EV1;
end

function u =get_utilities(thetas,data,states,EV,beta)
% This function gets the utilities of each alternative given the
% paramenters and the value functions. 
state = states(data);
u_replace = -thetas(1)*ones(length(data),1) + beta*EV(1);
u_notreplace = -0.001*thetas(2)*(state) + beta*EV(data);
u = [u_replace,u_notreplace];
end

%% 
