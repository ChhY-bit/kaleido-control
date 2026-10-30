clear,clc,close all

%% 参数设置
Ac = [0,1;-1,-2];
Bc = [0;1];
C = [1,0];
D = 1;

Ts = 0.05;

N = 30;
Q = 2;
R = 1;
S = 0.1;

%% 初始化
sys = c2d(ss(Ac,Bc,C,D),Ts);
A = sys.A;
B = sys.B;

controller = ctr.mpc(A,B,C,D,N,Q,R,S);