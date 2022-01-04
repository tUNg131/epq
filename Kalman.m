function Kalman(data, var_ran_walk, i)
% data = (angular_rate, temperature, cos_alpha)
[row, ~] = size(data);

n = row;
I = eye(3);
H = zeros(1, 3, n); %(cos, T, 1)
% X_( angular rate; temp coefficient; b)

Z = data(1:end,1);
H(1:end,1,1:end) = 1; %cos_alpha = 1 (parralel to the rotation axis)
H(1:end,2,1:end) = data (1:end,2);
H(1:end,3,1:end) = 1;
datalog = zeros(n,3);

% Initialization
delta = 0.001;
Q = [ delta 0     0;
    0       delta 0;
    0       0     var_ran_walk];
R = 1;


X_(1,1) = 0;
X_(2,1) = 0.05;
X_(3,1) = -1.8;

P_         = [
    100 0 0;
    0 100 0;
    0 0 100
    ]; 
% forgetting factor ( a larger puts more weights on previous estimates and therefore incurs less fluctuation of Rk, and longer time delays to catch up with changes)
% E_d = 0;
% E_e = 0;
% alpha = 0.995;

for k = 1:n-1
    
    
    
    d = Z(k) - H(:,:,k) * X_; %Innovation
    
    
    
%     E_d = (E_d * (k-1) + d * transpose(d))/k; % Averging d overtime
    
    K = P_*transpose(H(:,:,k))/(H(:,:,k)*P_*transpose(H(:,:,k)) + R); %Kalman gain
    
    X = X_ + K * d; %Correction
    
    e = Z(k) - H(:,:,k) * X; % Residual
    
    alpha = 0.995;
    
    if e >= 0.03 || e<= -0.03
        alpha = 0.1;
    end
%     alpha = 5.001 ^ (-e/0.01);
%     E_e = (E_e * (k-1) + e * transpose(e))/k; % Averging e overtime
    
    P = (I - K * H(:,:,k)) * P_; % State Covariance Update
    
    Q = alpha * Q +(1-alpha) * (K * d * transpose(d) * transpose (K)); % State Covariance Update
    
%     Q = K * E_d * transpose(K); % State Covariance Update
% 
%     Q(3,3) =  var_ran_walk;
%     Q(1:2,3) = 0;
%     Q(3,1:2) = 0;
    
    %---------------- k = k+1--------------
    
    X_ =  X; %State prediction
    
    P_ =  P + Q; % Priori Covariance Matrix
    
%     alpha = 1.001 ^ (-d/100);
    
    R =  R + (1-alpha)*(e * transpose(e) + H(:,:,k+1) * P_ * transpose(H(:,:,k+1))); %Measurement Covariance update
    
%     R = E_e + H(:,:,k+1) * P_ * transpose(H(:,:,k+1)) ; %Measurement Covariance update
    
    %---------------Data log---------------
    
    datalog (k,1) = k * i;
    datalog (k,2) = X_(1,1);
    datalog (k,3) = sqrt(Q(3,3));
end

% Plot

hold on;
plot(datalog (1:n-1,1) , datalog (1:n-1,2) ,'r');
plot(datalog (1:n-1,1) , datalog (1:n-1,2) - 2 * datalog (1:n-1,3), '-.b');
plot(datalog (1:n-1,1) , datalog (1:n-1,2) + 2 * datalog (1:n-1,3), '-.b');
yline(0.004,'-.g');
yline(-0.004,'-.g');
scatter((1:n-1)*i,Z(1:n-1),5,'r','filled');
% yline(0,'-.');
hold off;

% ylabel('G');
% title('Gravitational acceleration in Concord College');
end
