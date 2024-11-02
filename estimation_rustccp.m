clear all
close all
clc



load('cleandata.mat');


% First select the amount of buses

data = data(:,1:67);
choices = choices(:,1:67);
% GENERATE TRANSITION PROBABILITIES

% Compute increments


increment = data(2:end,:) - data(1:end-1,:);

% replace negative increments with nan

increment(increment<0) = NaN;

% compute frequencies

total = tabulate(reshape(increment.',1,[]));

trans = total(:,2)./sum(total(:,2));

% compute standard errors

standarderrors=sqrt(trans.*(1-trans)/(sum(total(:,2))));

% Estimation CCP

% define some parameters

beta=0.999;
states = linspace(1,90,90)';

% load the data
[T,N] = size(data);
trans_matrix = get_trans_matrix(trans);
theta1 =[0,0];
thetareal = [12, 2];




ccps = get_ccps(data,choices);
realccp = get_realccps(thetareal,states,beta,trans_matrix);
logitccps = logit_ccps(data,choices,states);

% estimate the likelihood
[theta13,fval13,exti,output,grad,hessian13]= fminunc(@(thetas) loglike_ccp(thetas,beta,choices,data,realccp,trans_matrix), theta1);

theta13

% let's construct the likelihood

function out = loglike_ccp(thetas,beta,choices,data,ccps,trans_matrix)


[T,N] = size(data);

% loop over all periods
ft =zeros(T,1);

% get expected ccps:

expected = trans_matrix*log(ccps);

for t=1:T

    % drop nans
    datanew = data(t,:);
    datanew = datanew(~isnan(datanew));

    choicesnew = choices(t,:);
    choicesnew = choicesnew(~isnan(choicesnew));
    
    loglike  =  zeros(N,1);
    for i=1:length(datanew)
        % get the payoff difference using ccps so that continuation value
        % cancels out. 
        diff  = -thetas(1) +beta*(-thetas(1) -log(ccps(1))) - (-0.001*thetas(2)*(datanew(i))+beta*(-thetas(1) -expected(datanew(i))));
        ccp_norepalce = 1/(1+exp(diff));
        % construct the likelihood
        loglike(i) =  choicesnew(i)*log(1-ccp_norepalce)...
                      +(1-choicesnew(i))*log(ccp_norepalce);
        
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

function predictedccps = logit_ccps(data,choices,states)

% generate a column vector

data = reshape(data.',1,[]);
choices = reshape(choices.',1,[]);

% drop nans

data = data(~isnan(data))';
choices = choices(~isnan(choices))';

% include squared

data(:,2) = data(:,1).^2;
states(:,2) = states(:,1).^2;

% fit a logit model

mdl = fitglm(data,choices);

% predict

predictedccps = predict(mdl,states);

predictedccps(predictedccps<0) = 0.01;

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


function realccp = get_realccps(thetareal,states,beta,trans_matrix)

% get continuation value
EV=valuefunction(states,trans_matrix,beta,thetareal);

% First get utilities
u = get_utilities(thetareal,states,EV,beta);

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