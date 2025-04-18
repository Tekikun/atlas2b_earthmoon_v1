clc;clear;

addpath('/Users/woodywu/Desktop/Research/General_Matlab');groot_setup;

Re = 6378.137; %km
AU = 149597870.700; %km

% file = 'atlas_pmax10_qmax10.dat';
% file = 'superatlasv3_massmoon_emoon.out';
file = 'superatlasv3.out';

% Open the file
fid = fopen(file,'r');
% Skip the first two lines of text
fgetl(fid);  % "I WILL CALCULATE..."
fgetl(fid);  % "pla  kp:k    a(au)   ..."
dataCell = textscan(fid, ...
    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'MultipleDelimsAsOne', true, ...
    'Delimiter', ' ');% Use textscan to read the remaining lines as 19 floating-point values
fclose(fid);% Close the file

% Convert the cell array to a numeric matrix (rows × 19 cols)
atlas_data = cell2mat(dataCell);
% Define new variable names for the first 11 columns
newVarNames = {'pla', 'kp', 'k', 'a_Re', 'e', 'i', 'w', 'ln', 'R_avg', 'R_diff', 'width_Re'};

%
% find uniq a in Re
a_list  = atlas_data(:,4);
e_list  = atlas_data(:,5);
a_width = atlas_data(:,11);
a_min = atlas_data(:,4) - a_width/2;
a_max = atlas_data(:,4) + a_width/2;

[a_uniq,ia,ic] = unique(a_list);
e_uniq = unique(e_list);

for idx = 1:1:length(a_uniq)
    a_tmp = a_uniq(idx);
    a_list_mat(:,idx) = a_list(a_list==a_tmp);
    a_min_mat(:,idx)  = a_min(a_list==a_tmp);
    a_max_mat(:,idx)  = a_max(a_list==a_tmp);
    p_list(idx,:) = atlas_data(ia(idx),2); % p value from data
    q_list(idx,:) = atlas_data(ia(idx),3); % q value from data
end

% choose range
p_list_select = p_list(p_list<=10 & q_list<=10);
q_list_select = q_list(p_list<=10 & q_list<=10);
a_uniq_select = a_uniq(p_list<=10 & q_list<=10);
% atlas string creation
atlas_n_strings = [num2str(p_list_select),repmat(':',length(p_list_select),1),num2str(q_list_select),...
                        repmat(char(9790),length(p_list_select),1)];

%
fig=figure('Units','inches','Position',[2,2,15,8]); hold on;
set( fig, 'Color', 'w', 'units', 'inches');
set( fig, 'position', [0,5,21,8]); 

% atlas by polygon
e_uniq_list = [e_uniq;e_uniq(end:-1:1)];
for i = 1:1:length(a_uniq)
pgone{i} = polyshape([a_min_mat(:,i);a_max_mat(end:-1:1,i)],e_uniq_list);
pg{i} = plot(pgone{i});
pg{i}.LineStyle = 'none';
% pg{i}.FaceColor = [192,192,192]/255;
end

% string
labelpoints(a_uniq_select,0.9,atlas_n_strings,'rotation',90,'FontSize',10)
% IC condition position
% xline(IC(1),'Color','r','LineWidth',2);


xlabel('Semi-major axis (Earth radii)');
ylabel('Eccentricity')

ylim([0,1])
xlim([10,90])

