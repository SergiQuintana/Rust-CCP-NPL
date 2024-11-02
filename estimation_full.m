
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
choices(choices==0) =2;
trans_matrix = get_trans_matrix(trans);
theta1 =[0,0];


% estimate the likelihood
[theta13,fval13,exti,output,grad,hessian13]= fminunc(@(thetas) loglike_full(thetas,beta,trans_matrix,choices,data,states), theta1);

theta13

% let's construct the likelihood

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


function C = get_continuation(data,EV,trans_matrix)

expected = trans_matrix*EV;
C  = ones(length(data),1);


for i=1:length(data)
C(i) = expected(data(i));
end

end

function vjt = get_conditional(data,theta,EV,beta,trans_matrix)
%This function draws a gumbell and perfomrs the chocie

% First get utilities
u = get_utilities(theta,data,EV,beta);

C = get_continuation(data,EV,trans_matrix);

% get conditional value function
v1 = u(:,1) + beta*EV(1);
v2 = u(:,2) + beta*C;
vjt = [v1,v2];

end

function u =get_utilities(thetas,data,EV,beta)
% This function gets the utilities of each alternative given the
% paramenters and the value functions. 
u_replace = -thetas(1)*ones(length(data),1);
u_notreplace = -0.001*thetas(2)*(data');
u = [u_replace,u_notreplace];
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

