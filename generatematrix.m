function [A_2, A_3, A_4, A_5, A_6, A_7]=generatematrix(coupling_type)

%Copyright (c) 2016 by Gary Iacobucci

%[A_2, A_3, A_4, A_5, A_6, A_7]=generatematrix(2);
%coupling_type, 1 = positive, 2 = negative

syms r z e d n k %adapted notation from Chung and Kennedy 1996

%Generate independent transition probability matrices for N channel patches
V = [z 1-z;... %single channel
     1-r r];
Pi_2 = kron(V,V); %two channels
Pi_3 = kron(V,Pi_2); %three channels
Pi_4 = kron(Pi_3,V); %four channels
Pi_5 = kron(Pi_4,V); %five channels
Pi_6 = kron(Pi_5,V); %six channels
Pi_7 = kron(Pi_6,V); %seven channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate aggregated matrices
L_2 = [1 0 0 0;... %3x4
0 2 2 0;...
0 0 0 1];
R_2 = [1 0 0;... %4x3
0 1 0;...
0 1 0;...
0 0 1];

R_3 = [1 0 0 0;... %8x4
0 1 0 0;...
0 1 0 0;...
0 1 0 0;...
0 0 1 0;...
0 0 1 0;...
0 0 1 0;...
0 0 0 1];
L_3 = [1 0 0 0 0 0 0 0;... %4x8
0 3 3 3 0 0 0 0;...
0 0 0 0 3 3 3 0;...
0 0 0 0 0 0 0 1];

L_4 = zeros(5,16);
L_4(1,1) = 1;
L_4(2,2:5)=4;
L_4(3,6:11) = 5;
L_4(4,12:15) = 4;
L_4(5,16) = 1;
R_4 = zeros(16,5);
R_4(1,1) = 1;
R_4(2:5,2) = 1;
R_4(6:11,3) = 1;
R_4(12:15,4) = 1;
R_4(16,5) = 1;

L=5;
R_5 = zeros(2^L,L+1); %32x6
R_5(1,1) = 1;
R_5(2:6,2) = 1;
R_5(7:16,3) = 1;
R_5(17:26,4) = 1;
R_5(27:31,5) = 1;
R_5(32,6) = 1;
L_5 = zeros(L+1,2^L); %6x32
L_5(1,1) = 1;
L_5(2,2:6)=5;
L_5(3,7:16) = 10;
L_5(4,17:26) = 10;
L_5(5,27:31) = 5;
L_5(6,32) =1;

L = 6;
R_6 = zeros(2^L,L+1); %32x6
R_6(1,1) = 1;
R_6(2:7,2) = 1;
R_6(8:22,3) = 1;
R_6(23:42,4) = 1;
R_6(43:57,5) = 1;
R_6(58:63,6) = 1;
R_6(64,7) = 1;
L_6 = zeros(L+1,2^L); %6x32
L_6(1,1) = 1;
L_6(2,2:7)=6;
L_6(3,8:22) = 15;
L_6(4,23:42) = 20;
L_6(5,43:57) = 15;
L_6(6,58:63) =6;
L_6(7,64) = 1;

L = 7;
R_7 = zeros(2^L,L+1); %32x6
R_7(1,1) = 1;
R_7(2:8,2) = 1;
R_7(9:29,3) = 1;
R_7(30:64,4) = 1;
R_7(65:99,5) = 1;
R_7(100:120,6) = 1;
R_7(121:127,7) = 1;
R_7(128,8) = 1;
L_7 = zeros(L+1,2^L); %6x32
L_7(1,1) = 1;
L_7(2,2:8)=6;
L_7(3,9:29) = 15;
L_7(4,30:64) = 20;
L_7(5,65:99) = 15;
L_7(6,100:120) =6;
L_7(7,121:127) = 1;
L_7(8,128) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate coupling matrices
switch coupling_type
    case 1 %positive coupling
        Pc_2 = sym(zeros(4,4));
        Pc_2(1,1) = z;
        Pc_2(1,end) = (1-z);
        Pc_2(end,1) = (1-r);
        Pc_2(end,end) = r;
        Pc_2(2:end-1,1) = 0.5;
        Pc_2(2:end-1,end) = 0.5;

        Pc_3 = sym(zeros(8,8));
        Pc_3(1,1) = z;
        Pc_3(1,end) = (1-z);
        Pc_3(end,1) = (1-r);
        Pc_3(end,end) = r;
        Pc_3(2:end-1,1) = 0.5;
        Pc_3(2:end-1,end) = 0.5;

        Pc_4 = sym(zeros(16,16));
        Pc_4(1,1) = z;
        Pc_4(1,end) = (1-z);
        Pc_4(end,1) = (1-r);
        Pc_4(end,end) = r;
        Pc_4(2:end-1,1) = 0.5;
        Pc_4(2:end-1,end) = 0.5;

        Pc_5 = sym(zeros(32,32));
        Pc_5(1,1) = z;
        Pc_5(1,end) = (1-z);
        Pc_5(end,1) = (1-r);
        Pc_5(end,end) = r;
        Pc_5(2:end-1,1) = 0.5;
        Pc_5(2:end-1,end) = 0.5;

        Pc_6 = sym(zeros(64,64));
        Pc_6(1,1) = z;
        Pc_6(1,end) = (1-z);
        Pc_6(end,1) = (1-r);
        Pc_6(end,end) = r;
        Pc_6(2:end-1,1) = 0.5;
        Pc_6(2:end-1,end) = 0.5;

        Pc_7 = sym(zeros(128,128));
        Pc_7(1,1) = z;
        Pc_7(1,end) = (1-z);
        Pc_7(end,1) = (1-r);
        Pc_7(end,end) = r;
        Pc_7(2:end-1,1) = 0.5;
        Pc_7(2:end-1,end) = 0.5;

    case 2 %negative coupling
        i = size(Pi_2);
        Pc_2 = sym(zeros(i(1),i(2)));
        Pc_2(1,1) = e;
        Pc_2(1,2) = (1-e);
        Pc_2(2,1) = (1-n);
        Pc_2(2,2) = n;
        Pc_2(3:end,1) = (1-d);
        Pc_2(3:end,2) = d;

        i = size(Pi_3);
        Pc_3 = sym(zeros(i(1),i(2)));
        Pc_3(1,1) = e;
        Pc_3(1,2) = (1-e);
        Pc_3(2,1) = (1-n);
        Pc_3(2,2) = n;
        Pc_3(3:end,1) = (1-d);
        Pc_3(3:end,2) = d;

        i = size(Pi_4);
        Pc_4=sym(zeros(i(1),i(2)));
        Pc_4(1,1) = e;
        Pc_4(1,2) = (1-e);
        Pc_4(2,1) = (1-n);
        Pc_4(2,2) = n;
        Pc_4(3:end,1) = (1-d);
        Pc_4(3:end,2) = d;

        i = size(Pi_5);
        Pc_5=sym(zeros(i(1),i(2)));
        Pc_5(1,1) = e;
        Pc_5(1,2) = (1-e);
        Pc_5(2,1) = (1-n);
        Pc_5(2,2) = n;
        Pc_5(3:end,1) = (1-d);
        Pc_5(3:end,2) = d;

        i = size(Pi_6);
        Pc_6=sym(zeros(i(1),i(2)));
        Pc_6(1,1) = e;
        Pc_6(1,2) = (1-e);
        Pc_6(2,1) = (1-n);
        Pc_6(2,2) = n;
        Pc_6(3:end,1) = (1-d);
        Pc_6(3:end,2) = d;

        i = size(Pi_7);
        Pc_7=sym(zeros(i(1),i(2)));
        Pc_7(1,1) = e;
        Pc_7(1,2) = (1-e);
        Pc_7(2,1) = (1-n);
        Pc_7(2,2) = n;
        Pc_7(3:end,1) = (1-d);
        Pc_7(3:end,2) = d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate final matrices
A_2 = (1-k)*L_2*Pi_2*R_2 + k*L_2*Pc_2*R_2;
A_3 = (1-k)*L_3*Pi_3*R_3 + k*L_3*Pc_3*R_3;
A_4 = (1-k)*L_4*Pi_4*R_4 + k*L_4*Pc_4*R_4;
A_5 = (1-k)*L_5*Pi_5*R_5 + k*L_5*Pc_5*R_5;
A_6 = (1-k)*L_6*Pi_6*R_6 + k*L_6*Pc_6*R_6;
A_7 = (1-k)*L_7*Pi_7*R_7 + k*L_7*Pc_7*R_7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save output
save('A_3.mat','A_3');
save('A_4.mat','A_4');
save('A_5.mat','A_5');
save('A_6.mat','A_6');
save('A_7.mat','A_7');
