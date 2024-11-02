
% Code to simulate rust data

clear all

close all
clc


%% 
% Set Parameters
rng(1)
beta=0.999;
thetas = [12, 2];
trans = [0.29,0.7,0.01];

states = linspace(1,90,90)';
n = 100;
T = 60;
[data,choices] = simulate_data(states,thetas,trans,beta,n,T);

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





