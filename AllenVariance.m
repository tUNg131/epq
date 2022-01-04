
function [Result_data,var_ran_walk,i] = AllenVariance(data, frequency) 
[row,~] = size(data);
% row = 10^6;
% Data ( angular_rate temperature cos_alpha )
% frequency = 100;
x = floor(log10(row)); 

% x is the x-range on the loglog graph = log10(row/(frequency)) -log10(1/frequency). 

n = x/0.25 + 1; % number of calculating sums

SortedData = zeros(row,n,2);
AllenVar = zeros(2,n);

for k = 1:n
    i = floor(10^(0.25*(k-1))); %i is the number of data in the bin
    m = floor(row/i); %number of bins
    SortedData (1,k,1) = mean(data((1:i+1),1));
    for j = 2:m
        SortedData (j,k,1) = mean(data((j-1)*i+1 : j*i,1));
        SortedData (j-1,k,2) = (SortedData(j,k,1) - SortedData(j-1,k,1))^2;
    end
    AllenVar(1,k) = sqrt(sum(SortedData(1:m,k,2))/(2*(m-1)));
    AllenVar(2,k) = i/frequency;
end

% AllenVar(AllenVar == 0) = [];

%Plot 

loglog(AllenVar(2,1:end),AllenVar(1,1:end),'-or');
grid on;

%Result
[~,M] = min (AllenVar(1,1:end));
% M = 4;

Result_data = SortedData(1:end,M,1);


i = floor(10^(0.25*(M-1))); %i is the number of data in the bin in result
m = floor(row/i); %number of bins in result

    for j = 1:m
    Result_data(j,2) = mean(data((j-1)*i+1 : j*i,2));
    end
    
%delete all zeros

Result_data( ~any(Result_data,2), :  ) = [];  %rows
Result_data( :, ~any(Result_data,1) ) = [];  %columns 
var_ran_walk = AllenVar(1,M)^2;
end
 