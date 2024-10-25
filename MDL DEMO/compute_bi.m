function [bi] = compute_bi(a0,i,N)

Hindx = setdiff([1:N],[i]);
A_Hindx = a0(:,Hindx);
A_tilde_i = A_Hindx(:,1:N-2)-A_Hindx(:,N-1)*ones(1,N-2);
bi = A_Hindx(:,N-1)-a0(:,i);
bi = (eye(N-1) - A_tilde_i*(pinv(A_tilde_i'*A_tilde_i))*A_tilde_i')*bi;
bi = bi/norm(bi);