
% PDE model for simulating single cell gene expressions involved in
% cell-cell communications

for sim=1:100
N=100; % lattice size

%% cellular scale initilization
SC1=rand(N,N); % sender cells 1
SC2=rand(N,N); % sender cells 2
RC=rand(N,N);  % receiver cells

SC1(SC1>0.95)=1;SC1(SC1<0.95)=0;
SC2(SC2>0.95)=1;SC2(SC2<0.95)=0;
RC(RC>0.95)=1;RC(RC<0.95)=0;
%% molecular scale initilization
% L1=0.01*rand(N,N); % intial vales
% L2=0.01*rand(N,N);
% L3=0.01*rand(N,N);
% L4=0.01*rand(N,N);
% L5=0.01*rand(N,N);

R1=rand(N,N).*RC; % only expressed in receiver cells (RC); Receptor gene expression was assumed to be fixed but not evolved.
R2=rand(N,N).*RC;

% TF1=0.01*rand(N,N).*RC;
% TF2=0.01*rand(N,N).*RC;
% TF3=0.01*rand(N,N).*RC;
% 
% TG1=0.01*rand(N,N).*RC;
% TG2=0.01*rand(N,N).*RC;
% TG3=0.01*rand(N,N).*RC;
% TG4=0.01*rand(N,N).*RC;

D1=0.1*1e-0; D2=0.2*1e-0; D3=0.3*1e-0; D4=0.1*1e-0; D5=0.3*1e-0; % diffusion coefficient
r11=abs(normrnd(0.5,0.1));r21=abs(normrnd(0.4,0.1));r31=abs(normrnd(0.3,0.1));r41=abs(normrnd(0.2,0.1));r51=abs(normrnd(0.1,0.1)); % prodcution/release rate
r12=abs(normrnd(0.1,0.01));r22=abs(normrnd(0.2,0.1));r32=abs(normrnd(0.3,0.1));r42=abs(normrnd(0.4,0.1));r52=abs(normrnd(0.5,0.1));
d1=0.1; d2=0.1; d3=0.1; d4=0.1; d5=0.1; % degradation coefficient

% M1=rand(N,N)+0.1;M2=rand(N,N)+0.1;M3=rand(N,N)+0.1;M4=rand(N,N)+0.1;M5=rand(N,N)+0.1; % maximal concentration

%% interaction parameters in the network
b11=1;b21=0;b31=1;b41=0;b51=1;  % Li-R1
b12=0;b22=0;b32=1;b42=1;b52=1;  % Li-R2
alpha11=rand;alpha21=0; % Ri-TF1
alpha12=rand;alpha22=rand; % Ri-TF2
alpha13=0;alpha23=rand; % Ri-TF3
beta1=rand;beta2=rand;beta3=rand; % degradation rate of TFs
mu11=rand;mu12=rand;mu13=0;mu14=0; % TF1-TGi
mu21=rand;mu22=rand;mu23=0;mu24=0; % TF2-TGi
mu31=0;mu32=0;mu33=rand;mu34=0.1; % TF3-TGi
gama1=rand;gama2=rand;gama3=rand;gama4=rand; % degradation rate of TGs

%% PDE numerical simulation 

%% Steady-state of ligands (elliptic PDEs)
f1=r11*(SC1>0)+r12*(SC2>0);
f2=r21*(SC1>0)+r22*(SC2>0);
f3=r31*(SC1>0)+r32*(SC2>0);
f4=r41*(SC1>0)+r42*(SC2>0);
f5=r51*(SC1>0)+r52*(SC2>0);

L1=FDM_elliptic(D1,d1,f1);
L2=FDM_elliptic(D2,d2,f2);
L3=FDM_elliptic(D3,d3,f3);
L4=FDM_elliptic(D4,d4,f4);
L5=FDM_elliptic(D5,d5,f5);

%% Steady state of TFs and TGs
        TF1=1/beta1.*((b11.*L1.*R1+b21.*L2.*R1+b31.*L3.*R1+b41.*L4.*R1+b51.*L5.*R1).*alpha11+(b12.*L1.*R2+b22.*L2.*R2+b32.*L3.*R2+b42.*L4.*R2+b52.*L5.*R2).*alpha21);
        TF2=1/beta2.*((b11.*L1.*R1+b21.*L2.*R1+b31.*L3.*R1+b41.*L4.*R1+b51.*L5.*R1).*alpha12+(b12.*L1.*R2+b22.*L2.*R2+b32.*L3.*R2+b42.*L4.*R2+b52.*L5.*R2).*alpha22);
        TF3=1/beta3.*((b11.*L1.*R1+b21.*L2.*R1+b31.*L3.*R1+b41.*L4.*R1+b51.*L5.*R1).*alpha13+(b12.*L1.*R2+b22.*L2.*R2+b32.*L3.*R2+b42.*L4.*R2+b52.*L5.*R2).*alpha23);
       
        TG1=1/gama1*(TF1*mu11+TF2*mu21+TF3*mu31);
        TG2=1/gama2*(TF1*mu12+TF2*mu22+TF3*mu32);
        TG3=1/gama3*(TF1*mu13+TF2*mu23+TF3*mu33);
        TG4=1/gama4*(TF1*mu14+TF2*mu24+TF3*mu34);


%% Spatial coordinates of cells
SC1_num=sum(sum(SC1(:,:)==1));SC2_num=sum(sum(SC2(:,:)==1));RC_num=sum(sum(RC(:,:)==1));

[x,y]=find(SC1(:,:)==1);
SC1_xy=[x,y];
[x,y]=find(SC2(:,:)==1);
SC2_xy=[x,y];
[x,y]=find(RC(:,:)==1);
RC_xy=[x,y];
Cell_xy=[SC1_xy;SC2_xy;RC_xy];

Cell_lables=[ones(SC1_num,1);2*ones(SC2_num,1);3*ones(RC_num,1)];  % 1-SC1; 2- SC2; 3-RC

%% Expression matrix
SC1_num=sum(sum(SC1(:,:)==1));SC2_num=sum(sum(SC2(:,:)==1));RC_num=sum(sum(RC(:,:)==1));
L1_2D=L1(:,:);L2_2D=L2(:,:);L3_2D=L3(:,:);L4_2D=L4(:,:);L5_2D=L5(:,:);
R1_2D=R1(:,:);R2_2D=R2(:,:);TF1_2D=TF1(:,:);TF2_2D=TF2(:,:);TF3_2D=TF3(:,:);TG1_2D=TG1(:,:);TG2_2D=TG2(:,:);TG3_2D=TG3(:,:);
EM=zeros(14,SC1_num+SC2_num+RC_num);
EM(1,1:SC1_num)=L1_2D(sub2ind([100,100],SC1_xy(:,1),SC1_xy(:,2)));
EM(2,1:SC1_num)=L2_2D(sub2ind([100,100],SC1_xy(:,1),SC1_xy(:,2)));
EM(3,1:SC1_num)=L3_2D(sub2ind([100,100],SC1_xy(:,1),SC1_xy(:,2)));
EM(4,1:SC1_num)=L4_2D(sub2ind([100,100],SC1_xy(:,1),SC1_xy(:,2)));
EM(5,1:SC1_num)=L5_2D(sub2ind([100,100],SC1_xy(:,1),SC1_xy(:,2)));

EM(1,SC1_num+1:SC1_num+SC2_num)=L1_2D(sub2ind([100,100],SC2_xy(:,1),SC2_xy(:,2)));
EM(2,SC1_num+1:SC1_num+SC2_num)=L2_2D(sub2ind([100,100],SC2_xy(:,1),SC2_xy(:,2)));
EM(3,SC1_num+1:SC1_num+SC2_num)=L3_2D(sub2ind([100,100],SC2_xy(:,1),SC2_xy(:,2)));
EM(4,SC1_num+1:SC1_num+SC2_num)=L4_2D(sub2ind([100,100],SC2_xy(:,1),SC2_xy(:,2)));
EM(5,SC1_num+1:SC1_num+SC2_num)=L5_2D(sub2ind([100,100],SC2_xy(:,1),SC2_xy(:,2)));

EM(6,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=R1_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(7,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=R2_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(8,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TF1_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(9,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TF2_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(10,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TF3_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));

EM(11,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TG1_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(12,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TG2_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(13,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TG3_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));
EM(14,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TG3_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));

%% record the simulation results for each slide
xlswrite('Lable_Coordinates_100_slides.xlsx',[Cell_lables,Cell_xy],sim) ;
xlswrite('ExpressionMatrix_100_slides.xlsx',[Cell_lables';EM],sim) ;
end