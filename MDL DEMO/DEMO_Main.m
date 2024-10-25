clear all; close all; clc
%% Setting
load GSE19380 % Gene Expression Data (brain disease, molecular)
X=MK'; X(1:16,:)= []; % remove the first 16 references (pure samples)
[M L] = size(X); Kmin=2; Kmax=9; % ground-truth K=4
IsSimplex=0; % the data violates the full-additivity in (A2) (i.e., NO simplex structure), requiring preprocessing (cf. Remark 2)

% load Cuprite_1997_150_150_3 % Remote Sensing Data (mining site, hyperspectral)
% [C d O_index]= RASF(X,9,10/225); X(:,O_index)= []; % remove 10 outliers
% [M,L] = size(X); Kmin=3; Kmax=11; % ground-truth K=9
% IsSimplex=1; % the data satisfies the full-additivity in (A2) (i.e., WITH simplex structure), without requiring preprocessing

%% MDL
[K_est code_length time] = MDL(X,Kmin,Kmax,IsSimplex);

%% Plot
plot(Kmin:Kmax,code_length(Kmin:Kmax),'Linewidth',2,'Marker','.','MarkerSize',20)
grid on
xlabel('Number of Sources')
ylabel('Code Length')