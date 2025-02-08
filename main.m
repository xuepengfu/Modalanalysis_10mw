clear
clc
set(0, 'defaulttextinterpreter','latex')  
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% load A matrix
data = readyaml(fullfile('00_IEA-10.0-198-RWT.BD1.sum.yaml'));
init_loc = data.Init_Nodes_E1; %initial location of the blade

node_x0 = data.Init_Nodes_E1(:, 1:3);%node coordinates
node_r0 = data.Init_Nodes_E1(:, 4:end);

% A matrix
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

Bnew = [];
for j = 1:10    
Bnew(:,(j-1)*6+1:(j-1)*6+3) = Ball(:,1+3*j:3+3*j);
Bnew(:,(j-1)*6+4:(j-1)*6+6) = Ball(:,34+3*j:36+3*j);
end

minusMinv = Bnew(61:120,:);
minusMinvK = A(61:120,1:60);
minusMinvC = A(61:120,61:120);

M = inv(minusMinv);
MM = (M+M')/2;
K = -inv(minusMinv)*minusMinvK;
KK = (K+K')/2;
C = -inv(minusMinv)*minusMinvC;
CC = (C+C')/2;

coe = C./K;
coed = diag(coe);

%% modal shape
%complex
[VecA, ValueA] = eig(A);

eigValues = diag(ValueA);
imagAbsValues = abs(imag(eigValues));
[~, sortIndex] = sort(imagAbsValues);
sortedEigValues = eigValues(sortIndex);
sortedVecA = VecA(:,sortIndex);
positiveImagIndex = imag(sortedEigValues) > 0;
negativeImagIndex = imag(sortedEigValues) <= 0;
finalEigValues = [sortedEigValues(positiveImagIndex); sortedEigValues(negativeImagIndex)];
finalVecA = [sortedVecA(:,positiveImagIndex) sortedVecA(:,negativeImagIndex)];

ValueA1 = diag(finalEigValues);
VecA1 = finalVecA;

% a1 = VecA1' * A * VecA1
% a2 = VecA1'* VecA1

frequencies_complex = abs(imag(diag(ValueA1))) / (2*pi);
frequencies_complex_withoutdamping = abs(diag(ValueA1))/(2*pi);
dampingratio_complex = abs(real(diag(ValueA1)) ./ (frequencies_complex_withoutdamping*2*pi));

comedgemode = VecA1(2:6:2+6*9,2);% complex edgewise modal shape
comedgemode = [0;comedgemode];

comedgemodeamp = abs(comedgemode);

%comedgemodeflap = VecA1(1:6:1+6*9,2);
%comedgemodeflapnew = abs(comedgemodeflap);

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

%simedgemodeflap = VecKMnew(1:6:1+6*9,2);
%simedgemodeflapnew = simedgemodeflap * 1 / (simedgemodeedge(end));

VecKMnew' * VecKMnew

%modal mass stiffness and damping
Mr = VecKMnew' * MM * VecKMnew;
%Mr(abs(Mr) < 1E-5) = 0;
Kr = VecKMnew' * KK * VecKMnew;
%Kr(abs(Kr) < 1E-5) = 0;
Cr = VecKMnew' * CC * VecKMnew;
%% damping power estimation
% power estimation of the first edgewise mode
% two methods here:
% 1. F*v dt
% 2. modal analysis

% method 1 (undamped modal shape)
second_mode_frequency=frequencies_complex_withoutdamping(2);
T = 1 / second_mode_frequency;
% \dot{x}(t)
simedgemodeall = VecKMnew(:,2);
zloc = init_loc(:,3);

t = linspace(0,T,100);
A = 1;
x_t = simedgemodeall*A/max(simedgemodeall(56)) * sin(2 * pi * second_mode_frequency * t);

v_t = diff(x_t, 1, 2) / (t(2)-t(1));
v_t = [v_t, v_t(:,end)];

f_t = zeros(size(v_t));
for i = 1:size(v_t, 2)
    f_t(:, i) = CC * v_t(:, i);  
end

P_t = zeros(size(v_t));
for i = 1:size(v_t, 2)
    P_t(:, i) = f_t(:, i) .* v_t(:, i); 
end
P_t = [zeros(6,size(P_t,2));P_t];

Ptnew = cell(1, 6); 
for i = 1:6
    Ptnew{i} = P_t(i:6:i+6*10, :); 
end

totalpowernew = cell(1, 6); 
for i = 1:6
    totalpowernew{i} = sum(Ptnew{i}, 1); 
end

%dissipation power of all DoFs
totalpowernewnew = cell(1, 6); 
for i = 1:6
    totalpowernewnew{i} = trapz(t, totalpowernew{i})/T; 
end

sum([totalpowernewnew{:}]);


total_power_per_time = sum(P_t, 1);  

total_energy = trapz(t, total_power_per_time);  
total_power = total_energy/(T)


%method1 (complex modal shape)

edge_mode_com = abs((VecA1(1:60, 2))'); 
edge_mode_com = edge_mode_com * 1 / edge_mode_com(56);

edge_mode_phase = -1 * angle(VecA1(1:60, 2)); %-1 for meets the python output

T = 1/frequencies_complex(2);
t = linspace(0,T,100);

x_tt = [];
for i = 1 :60
    for j = 1 :length(t)
x_tt(i,j) = edge_mode_com(i) * cos(2 * pi * frequencies_complex(2) * t(j) + edge_mode_phase(i));
    end
end

v_t = diff(x_tt, 1, 2) / (t(2)-t(1));
v_t = [v_t, v_t(:,end)];

f_t = zeros(size(v_t));
for i = 1:size(v_t, 2)
    f_t(:, i) = CC * v_t(:, i);  
end

P_t = zeros(size(v_t));
for i = 1:size(v_t, 2)
    P_t(:, i) = f_t(:, i) .* v_t(:, i); 
end
P_t = [zeros(6,size(P_t,2));P_t];

Ptnew = cell(1, 6); 
for i = 1:6
    Ptnew{i} = P_t(i:6:i+6*10, :); 
end

totalpowernew = cell(1, 6); 
for i = 1:6
    totalpowernew{i} = sum(Ptnew{i}, 1); 
end

%dissipation power of all DoFs
totalpowernewnew1 = cell(1, 6); 
for i = 1:6
    totalpowernewnew1{i} = trapz(t, totalpowernew{i})/T; 
end

sum([totalpowernewnew1{:}]);


total_power_per_time = sum(P_t, 1);  

total_energy = trapz(t, total_power_per_time);  
total_power = total_energy/(T)


% method2
Power_modal = (dampingratio_complex(2)*2*pi*(A/max(simedgemodeall(56)))^2 ...
    *1*(2*pi*frequencies_complex_withoutdamping(2))^2)/(T)




