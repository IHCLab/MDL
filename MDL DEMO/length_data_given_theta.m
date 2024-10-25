function [length_data] = length_data_given_theta( X, alpha_est, sigma_square_est, A_est )
[M L] = size(X);
[M N] = size(A_est);
how_many_hundreds= 100;
NoT = 100*how_many_hundreds; % NoT: number of trials
s = dirichlet_rndm(alpha_est,NoT,123)'; % this will give column sum = 1
z = A_est*s; % each column of z corresponds to a trial
number_of_inf = 0;
length_data = 0;
for j = 1:L
    % ==================== combination of for loop and matrix form ====================
    for ii=1:how_many_hundreds
        zz = z(:,1+100*(ii-1):100*(ii)) - X(:,j)*ones(1,100);
        yy = sum( zz.^2 )/(-2*sigma_square_est);
        log_fn_vec = yy - ( 0.5*M*log(2*pi*sigma_square_est) );
        p_jj(ii) = mean( exp(log_fn_vec) );
    end
    p_j = mean(p_jj);
    % ==================== combination of for loop and matrix form ====================
    if or(- log(p_j) == Inf, log(p_j) == Inf)
        number_of_inf = number_of_inf +1;
    else
        length_data = length_data- log(p_j);
    end
end

% use the following compensation strategy if too many data points have Inf or -Inf
if ((number_of_inf/L)>0.095), length_data= length_data*L/(L-number_of_inf); end