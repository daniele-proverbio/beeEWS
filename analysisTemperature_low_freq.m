%% Analysis of temperature data

clc; clear; close all;

%% 1 -- From Udine Dataset, load data
% Provide personal directory where data are saved

% Internal Temperature
data_DatasetInternalInfested = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_internal_infested');
data_DatasetInternalUninfested = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_internal_uninfested');
data_DatasetExternal = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_external');

days = datetime(data_DatasetInternalInfested(:,1),'ConvertFrom','excel') - calyears(1); % We are in 2018 actualyìly: 2019 in the original dataset is a typo

% Colony strength (rounded to integer, as fraction of bees may result from
% formulae but have little sense biologically)
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

% Days of sampling for colony strength
days_CS_sampling_infested = datetime(['16/08/2018'; '05/09/2018'; '01/10/2018'; '09/11/2018'; '10/12/2018'], 'InputFormat', 'dd/MM/yyyy');
days_CS_sampling_uninfested = datetime(['11/08/2018'; '06/09/2018'; '01/10/2018'; '09/11/2018'; '10/12/2018'], 'InputFormat', 'dd/MM/yyyy');


%% 2 -- Visualization
colors = lines(5);
markers = {'o', 'p', 'd', 's', '^'};

figure() % Infested
hold on
yyaxis left
for i = 2 : 6
    plot(days, data_DatasetInternalInfested(:,i), linewidth=1)
end 
plot(days, data_DatasetExternal(:,2), linewidth=1)
title("Infested apiaries")
ylabel("T (\circC)", fontsize=15)
colororder('default')
yyaxis right
for i = 1 : 5
    scatter(days_CS_sampling_infested, CS_infested(:,i), 'filled','Color',colors(i,:),'Marker',markers(i), 'SizeData',50)
end 
ylabel('Colony Strength', fontsize=20)
colororder('default')
legend(["5N", "8R", "13R", "14R", "CR/12", "External", "5N", "8R", "13R", "14R", "CR/12"], fontsize=14)


figure() % Uninfested
hold on
yyaxis left
for i = 2 : 6
    plot(days, data_DatasetInternalUninfested(:,i), linewidth=1)
end 
plot(days, data_DatasetExternal(:,2), linewidth=1)
legend(["6P", "7C", "10Bc", "VF", "VN", "External"], fontsize=14)
title("Uninfested apiaries")
ylabel("T (\circC)", fontsize=15)
colororder('default')
yyaxis right
for i = 1 : 5
    scatter(days_CS_sampling_uninfested, CS_uninfested(:,i), 'filled','Color',colors(i,:),'Marker',markers(i), 'SizeData',50)
end 
ylabel('Colony Strength', fontsize=20)
colororder('default')
legend(["6P", "7C", "10Bc", "VF", "VN","External","6P", "7C", "10Bc", "VF", "VN"], fontsize=14)

% From this visualization, we immediately see that temperature and colony
% strength do not really match with one another...


%% 3a -- Naive detreding
% Not meaningful: (a) Fluctuations on T:ext are greater than those of T_int
% -> a simple differencing would just give the (adjusted) inverse of T_ext;
% (b) T_int seems to be much more stationary

detrended_DatasetInternalInfested = data_DatasetInternalInfested(:,2:end) - data_DatasetExternal(:,2);

figure()
hold on
for i = 1 : 5
    plot(days, detrended_DatasetInternalInfested(:,i), linewidth=1)
end 
plot(days, data_DatasetExternal(:,2), linewidth=1)
legend(["5N", "8R", "13R", "14R", "CR/12", "External"])
title("Infested apiaries")
ylabel("T (\circC)")

% Test for linear relationship: in fact: internal temperature is minimally
% influenced by the external one. Only when the colony is already dieing (infested colonies),
% the temperature becomes the same as the external one (plus some offset
% associated with coibentation of the apiary)

figure()
hold on
for i = 2 : 6
    scatter(data_DatasetExternal(:,2), data_DatasetInternalInfested(:,i),'filled')
end 
legend(["5N", "8R", "13R", "14R", "CR/12"])
title("Infested apiaries")
ylabel("T_{ext} (\circC)")
xlabel("T_{int} (\circC)")


figure()
hold on
for i = 2 : 6
    scatter(data_DatasetExternal(:,2), data_DatasetInternalUninfested(:,i),'filled')
end 
legend(["6P", "7C", "10Bc", "VF", "VN"])
title("Infested apiaries")
ylabel("T_{ext} (\circC)")
xlabel("T_{int} (\circC)")


%% 3b -- Better to use a Granger causality test first, because it seems that T_int and T_ext do not really correlate (and T_ext does not cause T_int)

% First step: change point detection on the infested apiaries (do I
% identify when the colony start dieing?)

% Second step: segment data according to being before or after the changepoint
% Only for colonies that collapse; for uninfested, use the whole dataset


% On the infested
titles = ["5N", "8R", "13R", "14R", "CR/12"];
idx_changepoint = zeros(5,1);

h = zeros(5,2);  % first columns: before changepoint; second columns: after chaingepoint
pvalue = zeros(5,2); % as for h

h_u = zeros(5,1);  % for the uninfested (no changepoint)
pvalue_u = zeros(5,1); % as for h_u

figure()
for i = 1:5
    A = data_DatasetInternalInfested(:,i+1);
    A = A(~isnan(A));
    B = data_DatasetExternal(:,2);
    B = B(~isnan(A));
    days_changepoint = days(~isnan(A));

    % First step
    [idx_collapse, r] = findchangepts(A,MaxNumChanges=1);
    disp(days_changepoint(idx_collapse)-3);
    idx_changepoint(i) = idx_collapse - 3;

    subplot(5,1,i)
    plot(days_changepoint, A)
    xline(idx_collapse-3)
    ylabel("T(\circC)")
    title(titles(i))

    % Second step
    % Segment
    if idx_changepoint(i) > 0
        infestedBeforeChangepoint = A(1:idx_changepoint(i));
        infestedAfterChangepoint = A(idx_changepoint(i):end);
        externalBeforeChangepoint = B(1:idx_changepoint(i));
        externalAfterChangepoint = B(idx_changepoint(i):end);
    end

    % Make data stationary
    infestedBeforeChangepoint = price2ret(infestedBeforeChangepoint);
    infestedAfterChangepoint = price2ret(infestedAfterChangepoint);
    externalBeforeChangepoint = price2ret(externalBeforeChangepoint);
    externalAfterChangepoint = price2ret(externalAfterChangepoint);
    [h(i,1), pvalue(i,1)] = gctest(externalBeforeChangepoint, infestedBeforeChangepoint);
    [h(i,2), pvalue(i,2)] = gctest(externalAfterChangepoint, infestedAfterChangepoint);    
end
% Except for the final part of the collapsed apiary CR/12 (where
% temperature sets at the external temperature, plus a small offset, as the
% whole colony had died), there is no causality between external and
% internal temperature -> I can use temperatures as they are

% On the uninfested (whole time series)
for i = 1:5
    A = data_DatasetInternalUninfested(:,i+1);
    A = A(~isnan(A));
    B = data_DatasetExternal(:,2);
    B = B(~isnan(A));

    uninfested_stat = price2ret(A);
    external_stat = price2ret(B);

    [h_u(i), pvalue_u(i)] = gctest(external_stat, uninfested_stat );
end
% Uninfested, but also here with subdivision around changepoint
titles_u = ["6P", "7C", "10Bc", "VF", "VN"];
idx_changepoint_u = zeros(5,1);

h_u2 = zeros(5,2);  % first columns: before changepoint; second columns: after chaingepoint
pvalue_u2 = zeros(5,2); % as for h

figure()
for i = 1:5
    A = data_DatasetInternalUninfested(:,i+1);
    A = A(~isnan(A));
    B = data_DatasetExternal(:,2);
    B = B(~isnan(A));
    days_changepoint = days(~isnan(A));

    % First step
    [idx_collapse, r] = findchangepts(A,MaxNumChanges=1);
    disp(days_changepoint(idx_collapse)-3);
    idx_changepoint_u(i) = idx_collapse - 3;

    subplot(5,1,i)
    plot(days_changepoint, A)
    xline(idx_collapse-3)
    title(titles_u(i))
    ylabel("T(\circC)")

    % Second step
    % Segment
    if idx_changepoint_u(i) > 0
        uninfestedBeforeChangepoint = A(1:idx_changepoint_u(i));
        uninfestedAfterChangepoint = A(idx_changepoint_u(i):end);
        externalBeforeChangepoint = B(1:idx_changepoint_u(i));
        externalAfterChangepoint = B(idx_changepoint_u(i):end);
    end

    % Make data stationary
    uninfestedBeforeChangepoint = price2ret(uninfestedBeforeChangepoint);
    uninfestedAfterChangepoint = price2ret(uninfestedAfterChangepoint);
    externalBeforeChangepoint = price2ret(externalBeforeChangepoint);
    externalAfterChangepoint = price2ret(externalAfterChangepoint);
    [h_u2(i,1), pvalue_u2(i,1)] = gctest(externalBeforeChangepoint, uninfestedBeforeChangepoint);
    [h_u2(i,2), pvalue_u2(i,2)] = gctest(externalAfterChangepoint, uninfestedAfterChangepoint);    
end
% Here, I have three apiaries whose temperature depends on the external
% one: "6P", "10Bc", "VN". Apparently, they are less insulated and receive
% the influence of the external environment. Hence, I need to detrend these
% time series appropriately (see below)

%% 4 -- Detrend by external temperature
% Fit the apiaries that depend on the external environment with a plynomial fit: find a
% relationship to enable the detrending of the coarse and finely-sampled data
y1 = data_DatasetInternalUninfested(:,2);
y1 = y1(~isnan(y1));
y2 = data_DatasetInternalUninfested(:,4);
y2 = y2(~isnan(y2));
y3= data_DatasetInternalUninfested(:,6);
y3 = y3(~isnan(y3));
y4= data_DatasetInternalInfested(:,6);
y4 = y4(~isnan(y4));
x = data_DatasetExternal(:,2);
x = x(~isnan(y1));

[fitresult_6P, gof_6P] = fitTemp(x, y1);
[fitresult_10BC, gof_10BC] = fitTemp(x, y2);
[fitresult_VN, gof_VN] = fitTemp(x, y3);
[fitresult_CR, gof_CR] = fitTemp(x, y4);

% However, this polynomial detrending is not very good (R2 about 0.3);
% better to use a nonparametric detrending method: LOESS over expanding
% windows. It is carried out in a dedicated EWS_extract.m file