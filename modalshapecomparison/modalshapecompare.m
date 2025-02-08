clear
clc
set(0, 'defaulttextinterpreter','latex')  
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
%% load A matrix 
%revised to find the line with A then read A matrix
% fid = fopen('00_IEA-10.0-198-RWT.1.lin');
% lines = textscan(fid, '%s', 'Delimiter', '\n');
% fclose(fid);
% targetLine = find(strncmp(lines{1}, 'A:', 2), 1);
% fid = fopen('00_IEA-10.0-198-RWT.1.lin');
% AArray = textscan(fid, repmat('%f', 1, 120), 120, 'Delimiter', '\t', 'HeaderLines', targetLine, 'ReturnOnError', false);
% fclose(fid);

data = readyaml(fullfile('00_IEA-10.0-198-RWT.BD1.sum.yaml'));
init_loc = data.Init_Nodes_E1; %initial location of the blade

node_x0 = data.Init_Nodes_E1(:, 1:3);
node_r0 = data.Init_Nodes_E1(:, 4:end);


fid = fopen('00_IEA-10.0-198-RWT.1.BD1.lin');
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
targetLine = find(strncmp(lines{1}, 'A:', 2), 1);
fid = fopen('00_IEA-10.0-198-RWT.1.BD1.lin');
AArray = textscan(fid, repmat('%f', 1, 120), 120, 'Delimiter', '\t', 'HeaderLines', targetLine, 'ReturnOnError', false);
fclose(fid);

A = cell2mat(AArray);

%%load B matrix 
%dXdx: A ; dXdu B
fid = fopen('00_IEA-10.0-198-RWT.1.BD1.lin');
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
targetLine = find(strncmp(lines{1}, 'dXdu:', 4), 1);
fid = fopen('00_IEA-10.0-198-RWT.1.BD1.lin');
BBrray = textscan(fid, repmat('%f', 1, 306), 120, 'Delimiter', '\t', 'HeaderLines', targetLine, 'ReturnOnError', false);
fclose(fid);

B = cell2mat(BBrray);

Ball = B(:,19:84);
%Ball = B(:,1:60);

Bnew = [];
for j = 1:10    
Bnew(:,(j-1)*6+1:(j-1)*6+3) = Ball(:,1+3*j:3+3*j);
Bnew(:,(j-1)*6+4:(j-1)*6+6) = Ball(:,34+3*j:36+3*j);
end

%Bnewnew = B(:,1:60);

minusMinv = Bnew(61:120,:);
minusMinvK = A(61:120,1:60);
minusMinvC = A(61:120,61:120);


M1 = readmatrix('m.txt');
minusMinv1 = readmatrix('minus_m_inv.txt');
K1 = readmatrix('k.txt');
minusMinvK1 = readmatrix('minus_m_inv_k 1.txt');
K2 = load('kpython.mat');
K2 = K2.k;

M = inv(minusMinv);
MM = (M+M')/2;
K = -inv(minusMinv)*minusMinvK;
KK = (K+K')/2;
C = -inv(minusMinv)*minusMinvC;
CC = (C+C')/2;

% MM-M1
% K - K1
% K1 - K2
% K1 - K1'
% minusMinv1 - minusMinv
% minusMinvK-minusMinvK1
% issymmetric(M)
% issymmetric(K)
% diff =norm(minusMinvK - minusMinvK')/norm(minusMinvK);
% diff =norm(C - C')/norm(C);
% 
% diff =norm(K - K')/norm(K);
% max(abs(C - C')) / max(abs(C))

coe = C./K;
coed = diag(coe);
%% modal shape
%complex
[VecA, ValueA] = eig(A);

ValueA*VecA- VecA*ValueA

eigValues = diag(ValueA);
imagAbsValues = abs(imag(eigValues));
[~, sortIndex] = sort(imagAbsValues);
sortedEigValues = eigValues(sortIndex);
sortedVecA = VecA(:,sortIndex);

%aaa = AMatriinlin*sortedVecA - sortedVecA*diag(sortedEigValues);
positiveImagIndex = imag(sortedEigValues) > 0;
negativeImagIndex = imag(sortedEigValues) <= 0;

finalEigValues = [sortedEigValues(positiveImagIndex); sortedEigValues(negativeImagIndex)];
finalVecA = [sortedVecA(:,positiveImagIndex) sortedVecA(:,negativeImagIndex)];

ValueA1 = diag(finalEigValues);
VecA1 = finalVecA;

frequencies_complex = abs(imag(diag(ValueA1))) / (2*pi);
frequencies_complex_withoutdamping = abs(diag(ValueA1))/(2*pi);
dampingratio_complex = abs(real(diag(ValueA1)) ./ (frequencies_complex_withoutdamping*2*pi));

comedgemode = VecA1(2:6:2+6*9,2);
comedgemode = [0;comedgemode];

comedgemodeamp = abs(comedgemode);

comedgemodeflap = VecA1(1:6:1+6*9,2);
comedgemodeflapnew = abs(comedgemodeflap);

comedgemodeedge = VecA1(2:6:2+6*9,2);



%frequencies_complex_withoutdamping(2)
%based on K M
[VecKM, ValueKM] = eig(KK,MM);

eigKMValues = diag(ValueKM);
[~, sortKMIndex] = sort(eigKMValues);
sortedKMEigValues = eigKMValues(sortKMIndex);
sortedKMfres = sqrt(sortedKMEigValues)/(2*pi);

VecKMnew = VecKM(:,sortKMIndex);
VecKMnew = -1* VecKMnew;

simedgemodeedge = VecKMnew(2:6:2+6*9,2);
simedgemodeedgenew = simedgemodeedge * 1 / (simedgemodeedge(end));

simedgemodeflap = VecKMnew(1:6:1+6*9,2);
simedgemodeflapnew = simedgemodeflap * 1 / (simedgemodeedge(end));

Realedge=ComplexModeToRealMode(comedgemodeedge);
Realedgenew = Realedge * 1 / (Realedge(end));

Realflap=ComplexModeToRealMode(comedgemodeflap);
Realflapnew = Realflap * 1 / (Realedge(end));

Edgedtu = readmatrix('DTUedge.csv');
Flapdtu = readmatrix('DTUFlap.csv');
%% 
h = figure 
set(h, 'position', [100 100 1300 350]);
tiledlayout(1,2);
nexttile
plot(init_loc(2:end,3),Realedgenew,'k','LineWidth',1.5);
hold on
plot(init_loc(2:end,3),simedgemodeedgenew,'r:','LineWidth',1.5);
plot(Edgedtu(:,1),Edgedtu(:,2),'b--','LineWidth',1.5);

legend('Complex Modal Shape','Undamped Modal Shape','Grinderslev et al. (2022)','Location','northwest')
title('Edgewise part of 1st edgewise mode')
xlabel('Spanwise location [m]')
ylabel('Edgewise [-]')
set(gca, 'Fontname', 'Times New Roman','linewidth',0.75,'fontsize',22);
nexttile
plot(init_loc(2:end,3),Realflapnew,'k','LineWidth',1.5);
hold on
plot(init_loc(2:end,3),simedgemodeflapnew,'r:','LineWidth',1.5);
plot(Flapdtu(:,1),Flapdtu(:,2),'b--','LineWidth',1.5);
legend('Complex Modal Shape','Undamped Modal Shape','Grinderslev et al. (2022)','Location','northwest')
title('Flapwise part of 1st edgewise mode')
xlabel('Spanwise location [m]')
ylabel('Flapwise [-]')
set(gca, 'Fontname', 'Times New Roman','linewidth',0.75,'fontsize',22);

exportgraphics(h,'modalshapecom.pdf')

% y = chirp(t,f0,t1,f1,method)
% 
% 
% fs = 100; 
% T = 20 * 60; 
% t = 0:1/fs:T; 
% 
% f0 = 1/60; 
% f1 = 1/4; 
% t1 = T; 
% 
% y = chirp(t, f0, t1, f1, 'linear'); 
% y = y / max(abs(y)); 
% 
% plot(t, y);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Chirp Signal');
