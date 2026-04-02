%% Analysis of temperature data

clc; clear; close all;

%% From Udine Dataset, load data
% Temperature
data = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\sonde_2F\all_temperature.xlsx');

data_table = readtable('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\sonde_2F\all_temperature.xlsx'); 
days = datetime(data_table.Datetime,'InputFormat', 'dd/MM/yy HH.mm.ss');

% Colony strength (rounded to integer, as fraction of bees may result from
% formulae but has little sense biologically
CS_infested = round([14000 20821.9 25350.6 24591.6 17305.2;
14623.4 11739.2 21429.1 18873.8 12194.6;
11283.8 0 11385 12548.8 5388.9;
8602 0 6072 13662 2530;
5060 0 0 12650 0
]); % recreated matrix with same hive ordering as for temperatures ("5N", "8R", "13R", "14R", "CR/12")

CS_uninfested = round([22972.4 23984.4 22567.6 25173.5 19961.7;
15180 14724.6 16217.3 14699.3 17204;
13839.1 10170.6 8930.9 10195.9 13383.7;
17710 18722 15180 15180 17710;
17710 15686 13915 13156 15433
]); % recreated matrix with same hive ordering as for temperatures ("6P", "7C", "10Bc", "VF", "VN")

days_CS_sampling_infested = datetime(['16/08/2018'; '05/09/2018'; '01/10/2018'; '09/11/2018'; '10/12/2018'], 'InputFormat', 'dd/MM/yyyy');
days_CS_sampling_uninfested = datetime(['11/08/2018'; '06/09/2018'; '01/10/2018'; '09/11/2018'; '10/12/2018'], 'InputFormat', 'dd/MM/yyyy');


%% Visualization
end_vis = 6359; % Stop on 5th November as for the other dataset

colors = lines(5);
markers = {'o', 'p', 'd', 's', '^'};

figure() % Infested
hold on
yyaxis left
for i = 2 : 6
    plot(days(1:end_vis), data(1:end_vis,i), linewidth=1)
end 
title("Infested apiaries")
ylabel("T (\circC)", fontsize=20)
colororder('default')
yyaxis right
for i = 1 : 5
    scatter(days_CS_sampling_infested, CS_infested(:,i), 'filled','Color',colors(i,:),'Marker',markers(i), 'SizeData',50)
end 
ylabel('Colony Strength', fontsize=20)
colororder('default')
legend(["5N", "8R", "13R", "14R", "CR/12", "5N", "8R", "13R", "14R", "CR/12"], fontsize=14)


figure() % Uninfested
hold on
yyaxis left
for i = 7 : 11
    plot(days(1:end_vis), data(1:end_vis,i), linewidth=1)
end 
title("Uninfested apiaries")
ylabel("T (\circC)", fontsize=15)
colororder('default')
yyaxis right
for i = 1 : 5
    scatter(days_CS_sampling_uninfested, CS_uninfested(:,i), 'filled','Color',colors(i,:),'Marker',markers(i), 'SizeData',50)
end 
ylabel('Colony Strength', fontsize=20)
colororder('default')
legend(["6P", "7C", "10Bc", "VF", "VN","6P", "7C", "10Bc", "VF", "VN"], fontsize=14)

% From this visualization, we immediately see that temperature and colony
% strength do not really match with one another...




% There are noteworthy daily fluctuations, that need to be suppressed