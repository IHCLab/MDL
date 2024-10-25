%=====================================================================
% Programmer: Chia-Hsiang Lin (Steven)
% E-mail: chiahsiang.steven.lin@gmail.com
% Web: http://m105.nthu.edu.tw/~s105064538/index.html
% Date: January 14, 2015
%======================================================================

function [s_est] = ST_FCLS(X,A) % Steven's fast FCLS algorithm
N = size(A,2);
[M L] = size(X);
d = mean(A,2);
U = A-d*ones(1,N);
R = U*U';
[eV D] = eig(R);
C = eV(:,M-N+2:end);
Xd = C'*(X-d*ones(1,L)); 
Ad = C'*(A-d*ones(1,N));

% estimate b_1 to b_N (method 2: closed-form)  &  inner product constants
for i = 1:N-1
    b(:,i) = compute_bi(Ad,i,N);
    const(i,1) = b(:,i)'*Ad(:,i+1);
end
b(:,N) = compute_bi(Ad,N,N);
const(N,1) = b(:,N)'*Ad(:,1);

% source estimation
s_est = zeros(N,L);
for i = 1:N
    for j = 1:L
        s_est(i,j) = ( const(i,1) - b(:,i)'*Xd(:,j) )/( const(i,1) - b(:,i)'*Ad(:,i) );
    end
end

% now sum-to-one is satisfied; next non-negativety is handled
for i = 1:L
    useful_indx = ones(N,1);
    neg_indx = ( s_est(:,i) < 0 );
    num_neg_indx = sum(neg_indx);
    while num_neg_indx > 0
        useful_indx( find( s_est(:,i) <= 0 ) ) = 0;
        B = Ad(:, find(useful_indx==1) ); pp = Ad*s_est(:,i);
        v = size(B,2); BB=B(:,1:v-1)-B(:,v)*ones(1,v-1); pppp=pp-B(:,v);
        pppp = BB*( (BB'*BB)^(-1) )*BB'*pppp;
        pp = pppp+B(:,v);
        s_est(:,i) = zeros(N,1);
        for j = 1:N
            if useful_indx(j)==1
            s_est(j,i) = ( const(j,1) - b(:,j)'*pp )/( const(j,1) - b(:,j)'*Ad(:,j) );
            end
        end
        neg_indx = ( s_est(:,i) < 0 );
        num_neg_indx = sum(neg_indx);
        B=[]; pp=[]; vvv=[]; BB=[]; pppp=[];
    end
end