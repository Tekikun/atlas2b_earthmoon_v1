clc;clear;

addpath('/Users/woodywu/Desktop/Research/General_Matlab');groot_setup;

Re = 6378.137; %km
AU = 149597870.700; %km

% % IBEX IC
% IC = [1.619894838141151E+05,9.054106177359252E-01,1.197285453375567E+01,...
%     1.023846861023071E+02,3.696703201854587E+01,2.412048198408102E+02];
% IC(1) = IC(1) / Re; % Re
% file = 'IBEX_atlas2bres.dat';
% TESS IC 
IC = [1.671088887083566E-03*AU,6.380932362872992E-01,2.735501851350083E+01,...
    2.295025011248003E+02,5.488495273462943E+01,1.800061571716802E+02];
IC(1) = IC(1) /Re; %Re
file = 'TESS_atlas2bres.dat';

atlas_table = readtable(file);
atlas_data = table2array(atlas_table);


% find p+q and q less than 5
atlas_n = atlas_data(atlas_data(:,2)<=9 & atlas_data(:,3)<=9,:);
% atlas_n = atlas_data(atlas_data(:,2)<=20 & atlas_data(:,3)<=20,:);
atlas_n_strings = [num2str(atlas_n(:,2)),repmat(':',length(atlas_n),1),num2str(atlas_n(:,3)),...
                        repmat(char(9790),length(atlas_n),1)];

fig=figure;hold on;
set( fig, 'Color', 'w', 'units', 'inches');
set( fig, 'position', [0,5,21,8]); 
% atlas
stem(atlas_data(:,4),atlas_data(:,9),'Color','k','Marker', 'none');
% string
labelpoints(atlas_n(:,4),2*atlas_n(:,9),atlas_n_strings,'rotation',90)
% IC condition position
xline(IC(1),'Color','r','LineWidth',2);

set(gca, 'YScale', 'log')
xlim([5,80]);
% xlim([35,45])
ylim([1e-10,0.02]);

