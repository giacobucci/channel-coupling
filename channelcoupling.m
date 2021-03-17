function [output]=channelcoupling(filename,coupling_type,method,deadtime)    
%[output]=channelcoupling('GI-20161005a.dwt',2,1,1);
%[output]=channelcoupling('GI-20160810e.dwt',2,1,3);

%Copyright (c) 2016 by Gary Iacobucci
tic
fid=fopen(filename); %accepts DWT file from QUB
states = [];

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    celldata = textscan(tline,'%f %f');
    matdata = cell2mat(celldata);
    % match fails for text lines, textscan returns empty cells
    switch method
        case 1
            dt = 0.025*deadtime; %msec
            sample_points = ceil(matdata(:,2)/dt);
            idl_record=zeros(sample_points,1); %preallocate 1 columns 
            idl_record(:,1) = matdata(:,1); %all sample points equal to corresponding state
            states = [states; idl_record]; %append matrix to states
        case 2
            states = [states; matdata];
    end
end
       
fclose(fid);
syms z r e d n k;

[A_2, A_3, A_4, A_5, A_6, A_7]=generatematrix(coupling_type);

states = states(:,1)+1;
%Calculate estimated Transition and emision matrices
[TransitionEst, EmissionEst]=hmmestimate(states,states); %uses Baum-Welch
%algorithm to estimate transition and emissionmatrices.
    
A_hat=TransitionEst;
%A_hat = [0.1 .4 .5;...
%         .7 0.05 .295;...
%         0.45 .4 0.15];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
L=size(A_hat,1)-1 % %max # of channels
    
F=0;
if L==2
    for i=1:L
        for j=1:L
            F=F+(A_2(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==3
    for i=1:L
        for j=1:L
            F=F+(A_3(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==4
    for i=1:L
        for j=1:L
            F=F+(A_4(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==5
    for i=1:L
        for j=1:L
            F=F+(A_5(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==6
    for i=1:L
        for j=1:L
            F=F+(A_6(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==7
    for i=1:L
        for j=1:L
            F=F+(A_7(i,j)-A_hat(i,j))^2;
        end
    end
elseif L==8
    for i=1:L
        for j=1:L
            F=F+(A_Eight(i,j)-A_hat(i,j))^2;
        end
    end
end
    
    
F=0.5*F;    
%initial guess for [kappa,rho,zeta, eta, delta nu]
theta_o=[0.5 0.5 0.5 0.5 0.5 0.5];
    
%caculate the symbolic expression for the gradient
dFdtheta=[diff(F,k) diff(F,r) diff(F,z) diff(F,e) diff(F,d) diff(F,n)];
    
%initialize parameters gradient descent
thetaold=theta_o;
thetanew=[0 0 0 0 0 0];
precision=0.000005;
mu=0.001; %scaling factor for the gradient
    
loop=1;
iteration=1;
while( (norm(thetaold-thetanew)/norm(thetaold))>precision)
    thetaold=thetanew;
    thetanew=double(thetaold-mu*subs(dFdtheta,{k,r,z,e,d,n},{thetaold(1,1),thetaold(1,2),thetaold(1,3),thetaold(1,4),thetaold(1,5),thetaold(1,6)}))
    moment=thetanew(1,1)*thetanew(1,2)/thetanew(1,3);
    iteration=iteration+1;
    if iteration >=20000
        loop=loop+1;
        precision=precision*10;
        iteration=1;
    end 
    if loop==4
        break
    end
end
    
%display(iteration);
    
    %re-caculate if any element of theta is outside [0,1]
loop=1;  %initialize variable to prevent infinite loops
    %try to get kappa, rho, zeta within [0,1]
while( thetanew(1,1)<0 || thetanew (1,1)>1 || thetanew(1,2)<0 || ...
        thetanew(1,2)>1 || thetanew(1,3)<0 || thetanew(1,3)>1 || ...
        thetanew(1,4)>1 || thetanew(1,4)<0 || thetanew(1,5)>1 || ...
        thetanew(1,5)<0 || thetanew(1,6)>1 || thetanew(1,6)<0 || ...
        isnan(thetanew(1,1)) || isnan(thetanew(1,2)) || isnan(thetanew(1,3)) || isnan(thetanew(1,4)) || isnan(thetanew(1,5)) || isnan(thetanew(1,6)))
    precision=precision*1.001; %relax precision for every loop, 1.0001
    mu=mu/1.1; %decrease gradient scaling factor for every loop, 1.01
    thetaold=theta_o; %re-initialize variables
    thetanew=[0 0 0 0 0 0]; 
    iteration=1;
    loop=loop+1; %increment loop
    while( (norm(thetaold-thetanew)/norm(thetaold))>precision)
        thetaold=thetanew
        thetanew=double(thetaold-mu*subs(dFdtheta,{k,r,z,e,d,n},{thetaold(1,1),thetaold(1,2),thetaold(1,3),thetaold(1,4),thetaold(1,5),thetaold(1,6)}));
        moment=thetanew(1,1)*thetanew(1,2)/thetanew(1,3);
        iteration=iteration+1
        loop
        if iteration==1000
            break
        end
    end
        %only allow max 700 loops
    if loop==200
        break
    end
end
    
output.k=thetanew(1,1);
output.r=thetanew(1,2);
output.z=thetanew(1,3);
output.e=thetanew(1,4);
output.d=thetanew(1,5);
output.n=thetanew(1,6);
output.L=L;
output.moment=thetanew(1,1)*thetanew(1,2)/thetanew(1,3);
output.mu=mu;
output.precision=precision;
output.iteration=iteration;

toc
end