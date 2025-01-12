% PDE model for simulating single cell gene expressions involved in
% cell-cell communications, simulate slide-seqV2 or 10X visium data,
% multiple cells have same location.
clear all
N=100; % lattice size

base_save_path = 'F:/yll/博士相关/2.课题相关/4.stMLnet-yanllu/2024-11-24 GR修改/simulation_data/multiple_cells_simu_data/';

for sim=1:100
%% creat save path
save_path = fullfile(base_save_path, ['slide_' num2str(sim)]);
status = mkdir(save_path);
%% cellular scale initilization
SC1=rand(N,N); % sender cells 1
SC2=rand(N,N); % sender cells 2
RC=rand(N,N);  % receiver cells

SC1(SC1>0.95)=1;SC1(SC1<0.95)=0;
SC2(SC2>0.95)=1;SC2(SC2<0.95)=0;
RC(RC>0.95)=1;RC(RC<0.95)=0;

%% Spatial coordinates of cells
SC1_num=sum(sum(SC1(:,:)==1));SC2_num=sum(sum(SC2(:,:)==1));RC_num=sum(sum(RC(:,:)==1));

[x,y]=find(SC1(:,:)==1);
SC1_xy=[x,y];
[x,y]=find(SC2(:,:)==1);
SC2_xy=[x,y];
[x,y]=find(RC(:,:)==1);
RC_xy=[x,y];
Cell_xy=[SC1_xy;SC2_xy;RC_xy];
Cell_xy_0 = Cell_xy;
Cell_lables=[ones(SC1_num,1);2*ones(SC2_num,1);3*ones(RC_num,1)];  % 1-SC1; 2- SC2; 3-RC

%% setting x<x<x+10，y<y<y+10 areas have the same spatial coordiantes
patch_size=5;
for k=1:length(Cell_xy)
     x = Cell_xy(k,1);
     y = Cell_xy(k,2);
     new_x = x+patch_size;
     new_y = y+patch_size;
     condition_x = (Cell_xy(:, 1) > x(1)) & (Cell_xy(:, 1) < new_x(1));
     condition_y = (Cell_xy(:, 2) > y(1)) & (Cell_xy(:, 2) < new_y(1));
     result_indices = find(condition_x & condition_y);

     if ~isempty(result_indices)
            fprintf('满足条件的点索引：%s\n', mat2str(result_indices));
     end

     for l=1:length(result_indices)
         ind = result_indices(l);
         Cell_xy(ind,1) = x;
         Cell_xy(ind,2) = y;
     end
end

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


%% PDE numerical simulation 
D1=0.1*1e-0; D2=0.2*1e-0; D3=0.3*1e-0; D4=0.1*1e-0; D5=0.3*1e-0; % diffusion coefficient
r11=abs(normrnd(0.5,0.1));r21=abs(normrnd(0.4,0.1));r31=abs(normrnd(0.3,0.1));r41=abs(normrnd(0.2,0.1));r51=abs(normrnd(0.1,0.1)); % prodcution/release rate
r12=abs(normrnd(0.1,0.01));r22=abs(normrnd(0.2,0.1));r32=abs(normrnd(0.3,0.1));r42=abs(normrnd(0.4,0.1));r52=abs(normrnd(0.5,0.1));
d1=0.1; d2=0.1; d3=0.1; d4=0.1; d5=0.1; % degradation coefficient

%% interaction parameters in the network
b11=1;b21=0;b31=1;b41=0;b51=1;  % Li-R1
b12=0;b22=0;b32=1;b42=1;b52=1;  % Li-R2 ，根据图片的结构b52=0
alpha11=rand;alpha21=0; % Ri-TF1
alpha12=rand;alpha22=rand; % Ri-TF2
alpha13=0;alpha23=rand; % Ri-TF3
beta1=rand;beta2=rand;beta3=rand; % degradation rate of TFs
mu11=rand;mu12=rand;mu13=0;mu14=0; % TF1-TGi
mu21=rand;mu22=rand;mu23=0;mu24=0; % TF2-TGi
mu31=0;mu32=0;mu33=rand;mu34=0.1; % TF3-TGi
gama1=rand;gama2=rand;gama3=rand;gama4=rand; % degradation rate of TGs
%gama1=0.1 + (1 - 0.1) * rand;gama2=0.1 + (1 - 0.1) * rand;gama3=0.1 + (1 - 0.1) * rand;gama4=0.1 + (1 - 0.1) * rand;
%gama直接使用rand生成时，会出现1e-2的数值，导致1/gamma很大，进而导致基因表达量很大。

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

%% 
Cell_lables=[ones(SC1_num,1);2*ones(SC2_num,1);3*ones(RC_num,1)];  % 1-SC1; 2- SC2; 3-RC
colors = [1 0 0; 0 1 0; 0 0 1];  % 1:红色, 2:绿色, 3:蓝色
point_colors = colors(Cell_lables, :);  % 为每个点分配颜色

figure;
set(gcf, 'PaperPosition', [0, 0, 20, 28]); % 设置导出图像大小为 8x6 英寸
set(gcf, 'PaperSize', [20, 28]); % 设置纸张大小

subplot(1, 2, 1);
x = Cell_xy(:,1);y = Cell_xy(:,2);
scatter(x, y, 20, point_colors, 'filled');  
title('New spatial distribution');
xlabel('X Axis');
ylabel('Y Axis');

subplot(1, 2, 2);
x = Cell_xy_0(:,1);y = Cell_xy_0(:,2);
scatter(x, y, 20, point_colors, 'filled');  
title('Orginal spatial distribution');
xlabel('X Axis');
ylabel('Y Axis');
saveas(gcf, fullfile(save_path, 'Spatial_distribution.pdf'));

%% Expression matrix
SC1_num=sum(sum(SC1(:,:)==1));SC2_num=sum(sum(SC2(:,:)==1));RC_num=sum(sum(RC(:,:)==1));
L1_2D=L1(:,:);L2_2D=L2(:,:);L3_2D=L3(:,:);L4_2D=L4(:,:);L5_2D=L5(:,:);
R1_2D=R1(:,:);R2_2D=R2(:,:);TF1_2D=TF1(:,:);TF2_2D=TF2(:,:);TF3_2D=TF3(:,:);
TG1_2D=TG1(:,:);TG2_2D=TG2(:,:);TG3_2D=TG3(:,:);TG4_2D=TG4(:,:);
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
EM(14,SC1_num+SC2_num+1:SC1_num+SC2_num+RC_num)=TG4_2D(sub2ind([100,100],RC_xy(:,1),RC_xy(:,2)));

%% record the simulation results for each slide
xlswrite('Lable_Coordinates_100_slides.xlsx',[Cell_lables,Cell_xy],sim) ;
xlswrite('ExpressionMatrix_100_slides.xlsx',[Cell_lables';EM],sim) ;

%% plot gene expression vs spatial location

x = Cell_xy(:,1);y = Cell_xy(:,2); 

figure,
set(gcf, 'PaperPosition', [0, 0, 35, 28]); % 设置导出图像大小为 8x6 英寸
set(gcf, 'PaperSize', [35, 28]); % 设置纸张大小

subplot(2, 2, 1);
values = EM(11,:); % TG1 expression vs spatial location
scatter(x, y, 20, values, 'filled');  
colorbar;
title('Spatial gene expression of TG1');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG1_expression.png'));

subplot(2, 2, 2);
values = EM(12,:); % TG2 expression vs spatial location
scatter(x, y, 20, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG2');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG2_expression.png'));

subplot(2, 2, 3);
values = EM(13,:); % TG3 expression vs spatial location
scatter(x, y, 20, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG3');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG3_expression.png'));

subplot(2, 2, 4);
values = EM(14,:); % TG4 expression vs spatial location
scatter(x, y, 20, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG4');
xlabel('X Axis');
ylabel('Y Axis');
saveas(gcf, fullfile(save_path, 'FourGenes_Spatial_expression.pdf'));

end



