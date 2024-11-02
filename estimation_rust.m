
clear all
close all
clc

% CHANGE PATH HERE!
cd 'C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\Rust Simulation\rust data'

load('cleandata.mat');
%
%bus indices 1:15, 16:19, 20:67, 68:104
% define some parameters

beta=0.999;
states = linspace(1,90,90)';


% estimate everything

generate_table(beta,choices,data,states)

function generate_table(beta,choices,data,states)

% this function generates table 1 for Joan's notes. 

index13=[1,67];
index4=[68,104];
index14=[1,104];

disp('Estimating groups 1-3')
[trans13,sttrans13,thetahat13,stdtheta13] = get_estimates(index13,beta,choices,data,states);
disp('Groups 1-3 done!')
disp('Estimating group 4')
[trans4,sttrans4,thetahat4,stdtheta4] = get_estimates(index4,beta,choices,data,states);
disp('Group 4 done!')
disp('Estimating groups 1-4')
[trans14,sttrans14,thetahat14,stdtheta14] = get_estimates(index14,beta,choices,data,states);
disp('Groups 1-4 done!')

% Prepare tables frquencies

varnames = ["Estimates","Groups 1-3","Group 4","Groups 1-4"];
thetas = ["Phi 1";"Std Error";"Phi 2";"Std Error";"Phi 3";"Std Error"];

group13 = get_output_trans(trans13,sttrans13);
group4 = get_output_trans(trans4,sttrans4);
group14 = get_output_trans(trans14,sttrans14);

outputphi = table(thetas,group13,group4,group14,'VariableNames',varnames)
saveTableToLatex(outputphi, 'Output/transitions.tex')
% Prepare tables thetas

varnames = ["Estimates","Groups 1-3","Group 4","Groups 1-4"];
thetas = ["ThetaR";"Std Error";"Theta1";"Std Error"];

group13 = get_output_vector(thetahat13,stdtheta13);
group4 = get_output_vector(thetahat4,stdtheta4);
group14 = get_output_vector(thetahat14,stdtheta14);

outputtheta = table(thetas,group13,group4,group14,'VariableNames',varnames)
saveTableToLatex(outputtheta, 'Output/thetas_full.tex')

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
        row = T{i, :};
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

function out = get_output_trans(trans,std)

out = zeros(6,1);
out(1) = trans(1);
out(2) = std(1);
out(3) = trans(2);
out(4) = std(2);
out(5) = trans(3);
out(6) = std(3);
out = round(out,3);
end

function out = get_output_vector(thetahat,std)

out = zeros(4,1);
out(1) = thetahat(1);
out(2) = std(1);
out(3) = thetahat(2);
out(4) = std(2);
out = round(out,3);
end

function [trans,sttrans,thetahat,stdtheta] = get_estimates(index,beta,choices,data,states)


% First select the amount of buses

data = data(:,index(1):index(2));
choices = choices(:,index(1):index(2));

% compute transition frequencies

[trans,sttrans] = get_frequencies(data);

[T,N] = size(data);
choices(choices==0) =2;
trans_matrix = get_trans_matrix(trans);
theta1 =[0,0];

% estimate the model

[thetahat,fval13,exti,output,grad,hessianhat]= fminunc(@(thetas) loglike_full(thetas,beta,trans_matrix,choices,data,states), theta1);

% get st errors
stdtheta=sqrt(diag(inv(hessianhat)));

end

% construct the frequencies:

function [trans,sttrans] = get_frequencies(data)

% Compute increments


increment = data(2:end,:) - data(1:end-1,:);

% replace negative increments with nan

increment(increment<0) = NaN;

% compute frequencies

total = tabulate(reshape(increment.',1,[]));

trans = total(:,2)./sum(total(:,2));

% compute standard errors

sttrans=sqrt(trans.*(1-trans)/(sum(total(:,2))));

end

% let's construct the likelihood

function out = loglike_full(thetas,beta,trans_matrix,choices,data,states)

% Find the fixed point
EV = valuefunction(states,trans_matrix,beta,thetas);

[T,N] = size(data);

% loop over all periods
ft =zeros(T,1);

for t=1:T

        
    % drop nans
    datanew = data(t,:);
    datanew = datanew(~isnan(datanew));

    choicesnew = choices(t,:);
    choicesnew = choicesnew(~isnan(choicesnew));

    % If datanew is empty go to next period:

    if ~isempty(datanew)
        % get the conditional value functions
        
        vjt = get_conditional(datanew,thetas,EV,beta,trans_matrix);
        
        % set non-replacement as the base category: 
        
        vjt  = vjt - vjt(:,2);
        
        % get the ccps    
        ccps = exp(vjt)./(exp(vjt(:,1))+exp(vjt(:,2)));
        
        % use the ccp of the corresponding choice
        loglike  =  zeros(N,1);
        for i=1:length(datanew)
        
            loglike(i) = log(ccps(i,choicesnew(i)));
            
        end
        ft(t) = sum(loglike);

    end
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


