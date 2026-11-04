clear,clc
%%
Ts = 0.1;
A = [0,1,0;0,0,1;-1,-1,-1];
W = 0.1*diag([1,1,1]);
alpha = 10;

%%
A = expm(A*Ts);
M = dlyap(A',W);
[xt_min,xt_max,Et] = ctr.terminalSet_linear(A,M,alpha);

%%
% fig = figure(1);
% ctr.plots.plotEllipsoid(M,alpha,[],fig);
% ctr.plots.plotPolytope(Et,xt_min,xt_max,[],fig);