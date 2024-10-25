%=====================================================================
% Programmer: Chia-Hsiang Lin (Steven)
% E-mail: chiahsiang.steven.lin@gmail.com
% Web: http://m105.nthu.edu.tw/~s105064538/index.html
% Date: January 12, 2018
% -------------------------------------------------------
% Reference:
% C.-H. Lin, C.-Y. Chi, L. Chen, D. J. Miller, and Y. Wang
% ``Detection of sources in non-negative blind source separation by minimum description length criterion,"
% accepted by IEEE Trans. Neural Networks and Learning Systems, 2017.
%======================================================================
% A generic MDL algorithm for estimating the number of sources in non-negative blind source separation
% [K_est, code_length, time] = MDL(X,Kmin,Kmax,IsSimplex)
%======================================================================
%  Input
%  X is M-by-L data matrix, where M is the number of samples and L is the data length.
%  Kmin is the minimum number of sources.
%  Kmax is the maximum number of sources.
%  IsSimplex=1 if the full-additivity in assumption (A2) is satisfied; otherwise, IsSimplex=0.
%----------------------------------------------------------------------
%  Output
%  K_est is the estimated number of sources.
%  code_length is the code length of data, i.e., MDL(K).
%  time is the computation time (in secs).
%========================================================================

function [K_est code_length time] = MDL(X,Kmin,Kmax,IsSimplex)
load seed
t1 = clock;
[M L] = size(X);

% pre-processing & SNR pre-estimation by Corollary 1 (assuming K=Kmax)
% ================== Corollary 1 ==================
Upre = X - mean(X,2)*ones(1,L);
[eVpre Dpre] = eig( (Upre*Upre')/L );
lambda_vector_pre = Dpre*ones(M,1);
sigma_square_Kmax= sum( lambda_vector_pre(1:M-Kmax+1) )/(M-Kmax+1);
% ================== Corollary 1 ==================
snraa= (norm(X,'fro'))^2;
snrbb= M*L*sigma_square_Kmax;
SNR_pre_est= 10*log10(snraa/snrbb);

if IsSimplex==0, % making simplex structure of (2) violated, leading to invalid Gaussian-Dirichlet density
    X = X./(ones(M,1)*sum(X)); % standardized data (on F = aff{e_1,...,e_M})
    % =========== do PCA to present X in M-1 dimensional F ===========
    ddd = mean(X,2);
    UUU = X-ddd*ones(1,L);
    [eeVV DDD] = eig(UUU*UUU');
    CCC = eeVV(:,2:end);
    Xd = CCC'*(X-ddd*ones(1,L));
    X = [];
    X = Xd; Xd = [];
    [M L ] = size(X);
    % =========== do PCA to present X in M-1 dimensional F ===========
end

if M>100, % making density g of (21) almost infinite, leading to the numerical error MDL=-INF
    % =========== do dimension reduction to weaken g ===========
    M= 100;
    ddd = mean(X,2);
    UUU = X-ddd*ones(1,L);
    [eeVV DDD] = eig(UUU*UUU');
    CCC = eeVV(:,end-M+1:end);
    Xd = CCC'*(X-ddd*ones(1,L));
    X = [];
    X = Xd; Xd = [];
    [M L ] = size(X);
    % =========== do dimension reduction to weaken g ===========
end

if SNR_pre_est>40, % making density g of (21) almost 0, leading to the numerical error MDL=+INF
    % =========== add white noise to enhance g ===========
    SNR= 20; OBS= X; [M,L]= size(X); seed= 1;
    varianc = sum(OBS(:).^2)/10^(SNR/10) /M/L ;
    Cn = diag( varianc*ones(M,1) );
    randn('seed',seed);
    n = sqrtm(Cn)*randn([M L]);
    X = OBS+n;
    % =========== add white noise to enhance g ===========
end

length_theta = zeros(1,Kmax);
length_data = zeros(1,Kmax);
length_total = zeros(1,Kmax);
iiddxxOFoutlier = zeros(Kmax,L);

for K = 1:Kmin-1
fprintf('K = %d, ',K)
length_total(K)= Inf;
fprintf('MDL = %2.4f.\n',length_total(K))
end

for K = Kmin:Kmax
fprintf('K = %d, ',K)
N=K;

% ML estimation for sigma-square
U = X - mean(X,2)*ones(1,L);
[eV D] = eig( (U*U')/L );
lambda_vector = D*ones(M,1);
sigma_square_est = sum( lambda_vector(1:M-N+1) )/(M-N+1);

% ML estimation for a1 to aK
[A_est] = ADMM(X,N,seed(1));

% ML estimation for alpha
S_est = ST_FCLS(X,A_est);
[alpha_est] = dirichlet_fit_newton(S_est');

% coding length computation
length_theta(K) = 0.5*( 1+ ( N*(M+1) ) )*log(L);
length_data(K) = length_data_given_theta( X, alpha_est, sigma_square_est, A_est );
length_total(K) = length_theta(K) + length_data(K);
fprintf('MDL = %2.4f.\n',length_total(K))
end
time = etime(clock,t1);
fprintf('Computational Time = %2.4f (seconds).\n',time)
[~,K_est] = min(length_total);
fprintf('Estimated Number of Sources = %d.\n',K_est)
code_length = length_total;