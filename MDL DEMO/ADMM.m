% [49] J. M. Bioucas-Dias, ¡§A variable splitting augmented Lagrangian approach to linear spectral unmixing,¡¨ 
% in Proc. IEEE WHISPERS, Grenoble, France, Aug. 26-28, 2009, pp. 1¡V4.

function [MM]= ADMM(X,N,seed)
[MM,~,~,~] = sisal(X,N,seed);