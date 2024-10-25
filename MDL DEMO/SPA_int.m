function [A_est] = SPA_int(X,L,N)
[M L ] = size(X);
d = mean(X,2);
U = X-d*ones(1,L);
[eV D] = eig(U*U');
C = eV(:,M-N+2:end);
Xd = C'*(X-d*ones(1,L));
alpha_tilde = SPA(Xd,L,N); % the identified purest pixels
A_est = C * alpha_tilde + d * ones(1,N);

%% subprogram 1
function [bi] = compute_bi(a0,i,N)
Hindx = setdiff([1:N],[i]);
A_Hindx = a0(:,Hindx);
A_tilde_i = A_Hindx(:,1:N-2)-A_Hindx(:,N-1)*ones(1,N-2);
bi = A_Hindx(:,N-1)-a0(:,i);
bi = (eye(N-1) - A_tilde_i*(pinv(A_tilde_i'*A_tilde_i))*A_tilde_i')*bi;
bi = bi/norm(bi);
return;

%% subprogram 2
function [alpha_tilde] = SPA(Xd,L,N)

% Reference:
% [1] W.-K. Ma, J. M. Bioucas-Dias, T.-H. Chan, N. Gillis, P. Gader, A. J. Plaza, A. Ambikapathi, and C.-Y. Chi, 
% ``A signal processing perspective on hyperspectral unmixing,би 
% IEEE Signal Process. Mag., vol. 31, no. 1, pp. 67бV81, 2014.
% 
% [2] S. Arora, R. Ge, Y. Halpern, D. Mimno, A. Moitra, D. Sontag, Y. Wu, and M. Zhu, 
% ``A practical algorithm for topic modeling with provable guarantees,би 
% arXiv preprint arXiv:1212.4777, 2012.
%======================================================================
% An implementation of successive projection algorithm (SPA) 
% [alpha_tilde] = SPA(Xd,L,N)
%======================================================================
%  Input
%  Xd is dimension-reduced (DR) data matrix.
%  L is the number of pixels.   
%  N is the number of endmembers.
%----------------------------------------------------------------------
%  Output
%  alpha_tilde is an (N-1)-by-N matrix whose columns are DR purest pixels.
%======================================================================
%======================================================================

%----------- Define default parameters------------------
con_tol = 1e-8; % the convergence tolence in SPA
num_SPA_itr = N; % number of iterations in post-processing of SPA
N_max = N; % max number of iterations

%------------------------ initialization of SPA ------------------------
A_set=[]; Xd_t = [Xd; ones(1,L)]; index = [];
[val ind] = max(sum( Xd_t.^2 ));
A_set = [A_set Xd_t(:,ind)];
index = [index ind];
for i=2:N
    XX = (eye(N_max) - A_set * pinv(A_set)) * Xd_t;
    [val ind] = max(sum( XX.^2 )); 
    A_set = [A_set Xd_t(:,ind)]; 
    index = [index ind]; 
end
alpha_tilde = Xd(:,index);

%------------------------ post-processing of SPA ------------------------
current_vol = det( alpha_tilde(:,1:N-1) - alpha_tilde(:,N)*ones(1,N-1) );
for jjj = 1:num_SPA_itr
    for i = 1:N
        b(:,i) = compute_bi(alpha_tilde,i,N);
        b(:,i) = -b(:,i);
        [const idx] = max(b(:,i)'*Xd);
        alpha_tilde(:,i) = Xd(:,idx);
    end
    new_vol = det( alpha_tilde(:,1:N-1) - alpha_tilde(:,N)*ones(1,N-1) );
    if (new_vol - current_vol)/current_vol  < con_tol
        break;
    end
end
return;