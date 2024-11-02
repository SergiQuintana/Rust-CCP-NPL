%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       Rust (1987) Replication       %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code replicates Rust (1987) paper for particular bus group.
% It estimates the replacement and maintenance costs for bus groups 1,2 and
% 3. Estimation of different bus group(s) requires certain adjustments,
% commented in the code. 

% The script uses two other functions:
% - valuefunction.m - performs value function iteration algorithm
% - loglik_rust.m - consructs log-likelihood function for the Nested Fixed
%   Point Algorithm

% The following code is the solution of Problem Set 2 for the course of
% APPLIED MICROECONOMETRICS. Part I: Structural Micro, Spring 2024.

%Part I prepares the cleandata that was shared for the PS, starting from 
% the raw dataset. PS solution starts in Part II.

%% PART 1: CLEAN THE DATA

%% First step - Read and Upload the dataset 
%% 

clc; clear all; close all; 


% Use 'dlmread' command to read the data

rt50 = dlmread('rt50.asc');

t8h203 = dlmread('t8h203.asc');

g870 = dlmread('g870.asc');

a452372 = dlmread('a452372.asc');

a452374 = dlmread('a452374.asc');

a530872 = dlmread('a530872.asc');

a530874 = dlmread('a530874.asc');

a530875 = dlmread('a530875.asc');


% Extract relevant information
% Separate information data from observation for cumulative mileage

info1=   g870(1:11,1:end); g1=   g870(12:end,1:end);
info2=   rt50(1:11,1:end); g2=   rt50(12:end,1:end);
info3= t8h203(1:11,1:end); g3= t8h203(12:end,1:end);
info4=a530875(1:11,1:end); g4=a530875(12:end,1:end);
info5=a530874(1:11,1:end); g5=a530874(12:end,1:end);
info6=a452374(1:11,1:end); g6=a452374(12:end,1:end);
info7=a530872(1:11,1:end); g7=a530872(12:end,1:end);
info8=a452372(1:11,1:end); g8=a452372(12:end,1:end);

info=[info1,info2,info3,info4,info5,info6,info7,info8];

% Add information on date of replacement
% Add information about dates when replacements occur

start=74;

for i=1:size(info,2)

     info(12,i)= (info(11,i) -(start+1))*12 + info(10,i) ;
    
 if  info(4,i)==0
     info(13,i)=0;
 else
     info(13,i)=(info(5,i) -(info(11,i)+1))*12 + info(4,i) + (12 - info(10,i) +1);
 end
 if  info(7,i)==0
     info(14,i)=0;
 else
     info(14,i)=(info(8,i) -(info(11,i)+1))*12 + info(7,i) + (12 - info(10,i) +1);
 end
end

% Generate variable for replacement 
% Generate a matrix of decisions - taking values equal to 1 if a replecement 
% occurred. Note: it takes value NaN if bus was not in the fleet in a given
% month.

it=zeros(126,size(info,2));
for i=1:size(info,2)
    
    it(1:info(12,i),i)=NaN;
    
    if info(13,i)>0
        it(info(12,i)+info(13,i),i)=1;
    else
    end
    
    if info(14,i)>0
        it(info(12,i)+info(14,i),i)=1;
    else
    end
     
end

% Generate variable for cumulative mileage for BUS
% Combine cumulative mileage observations taking into account censored
% data

x1=[zeros((info1(11,1) -(start+1))*12 + info1(10,1),size(info1,2)); g1];
x2=[zeros((info2(11,1) -(start+1))*12 + info2(10,1),size(info2,2)); g2];
x3=[zeros((info3(11,1) -(start+1))*12 + info3(10,1),size(info3,2)); g3];
x4=[zeros((info4(11,1) -(start+1))*12 + info4(10,1),size(info4,2)); g4];
x5=[zeros((info5(11,1) -(start+1))*12 + info5(10,1),size(info5,2)); g5];
x6=[zeros((info6(11,1) -(start+1))*12 + info6(10,1),size(info6,2)); g6];
x7=[zeros((info7(11,1) -(start+1))*12 + info7(10,1),size(info7,2)); g7];
x8=[zeros((info8(11,1) -(start+1))*12 + info8(10,1),size(info8,2)); g8];

% Drop the first two observations for group 2 to make it the same
%size

x2=x2(3:end,:);  

% Adjust it for x2 accordingly! Basically move everything two periods up

it(78:end-2,16:19) = it(80:end,16:19);
it(end-2:end,16:19) = 0;


% Combine data and assign NaN for missing observation
xt=[x1,x2,x3,x4,x5,x6,x7,x8];
xt(find(xt==0))=NaN;


% Generate variable for monthly mileage of a BUS
% Take first difference of cumulative mileage to generate monthly
% mileage

dxt=zeros(size(xt,1),size(xt,2));
dxt(1,:)=NaN;
for i=1:size(xt,2)
    
    dxt(2:end,i)= xt(2:end,i)-xt(1:end-1,i);

end

% Check if there are any peculiarities

if find(dxt<0)
    display('something is wrong')
else
end


% Generate variable for cumulative mileage of ENGINE
% Combine cumulative mileage observations taking into account whether a
% replecement happened

% Generate vectors that take value 1 at the time of the replacement and
% every month afterwards

control1=zeros(size(xt,1),size(xt,2));
control2=zeros(size(xt,1),size(xt,2));

for i=1:size(xt,2)
    
    if    info(13,i)==0
        control1(:,i)=0;    
        control2(:,i)=0;
        
    elseif info(13,i)~=0 && info(14,i)==0
        
        control1(info(12,i)+info(13,i)+1:end,i)=  1;  
        control2(:,i)=0;
    elseif info(13,i)~=0 && info(14,i)~=0
        control1(info(12,i) + info(13,i): info(12,i) + info(14,i)-1,i)=  1; 
        control2(info(12,i) + info(14,i): end,i)=  1; 
    else
    end
end

% Construct matrix with mileage per engine, substracting the mileage at the
% replacement any month after the replacement

for i=1:size(xt,2)
    
    Ext(:,i)=xt(:,i) - control1(:,i).*info(6,i)  - control2(:,i).*info(9,i)    ;
    
end

% Replace negative values with original values as Rust (1987) did. 
% Take the mileage of replacement for true and disregarding the date when
% replacement happened - correct for all cases when Ext is negative.

EExt=Ext;

% Consider all the cases: mileage per engine can be negative in the month
% of the replacement, as well as one or two months afterwards

for i=1:size(xt,2)
        if info(14,i)==0
            EExt(find(Ext(:,i)<0),i)=xt(find(Ext(:,i)<0),i);
        else
            if find(Ext(info(12,i) + info(13,i),i)<0) >0
               EExt(info(12,i) + info(13,i),i)=xt(info(12,i) + info(13,i),i) ;
            else
            end
            
            if find(Ext(info(12,i) + info(13,i)+1,i)<0) >0
               EExt(info(12,i) + info(13,i)+1,i)=xt(info(12,i) + info(13,i)+1,i) ;
            else 
            end
            
            if find(Ext(info(12,i) + info(14,i),i)<0) >0
               EExt(info(12,i) + info(14,i),i)=xt(info(12,i) + info(14,i),i) - info(6,i);
            else
            end  
            if find(Ext(info(12,i) + info(14,i)+1,i)<0) >0
               EExt(info(12,i) + info(14,i)+1,i)=xt(info(12,i) + info(14,i)+1,i) - info(6,i);
            else  
            end
            if find(Ext(info(12,i) + info(14,i)+2,i)<0) >0
               EExt(info(12,i) + info(14,i)+2,i)=xt(info(12,i) + info(14,i)+2,i) - info(6,i);
            else  
            end
            
        end
end



% Generate a new dummy for replecement taking into account last changes
% Correct the decision matrix taking into account previous changes in the
% matrix with mileage per engine. Note that in cases with negative mileage, 
% we assume that the mileage at the replacement was correct and the date of
% replacement was wrong. 

EEit=zeros(size(EExt,1),size(EExt,2));

for i=1:size(EEit,2)
    for t=1:size(EEit,1)-1
        
        if EExt(t+1,i)<EExt(t,i)
        EEit(t,i)=1;
        else
        end
    end
end

for i=1:size(info,2)
    EEit(1:info(12,i),i)=NaN;
end

gridnumber=90; % number of states
upper=450000;
interval=upper/gridnumber;
state=[interval:interval:upper;1:1:gridnumber]';

for i=1:size(xt,2)
    for t=1:size(xt,1)
       
        differ=state(:,1)-EExt(t,i).*ones(size(state,1),1);
        
        differ(find(differ<0))=10000000;
        
        if differ(1)>0
             
             minimum=find(differ==min(differ));
             position(t,i)=state(minimum,2);
        else
             position(t,i)=NaN;
        end
        
    end   
end

% SERGI: Adjust it for gropu 2 accordingly! Basically move everything two periods up

EEit(78:end-2,16:19) = EEit(80:end,16:19);
EEit(end-2:end,16:19) = 0;
choices = EEit;
data = position;


save('cleandata',"data","choices")