%% Extract statistical indicators that may be used as EWS
% Here, from daily temperatures. 
% NB: I don't have a ground truth on when extintion happens: I can only check for false positives
% and negatives

clc; clear; close all

%% 1 -- Daily-sampled average temperatures
% Load data
% Update the folder where you save the dataset

data_DatasetInternalInfested = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_internal_infested');
data_DatasetInternalUninfested = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_internal_uninfested');
data_DatasetExternal = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\data_temperature\Dataset_external');

days = datetime(data_DatasetInternalInfested(:,1),'ConvertFrom','excel') - calyears(1); % We are in 2018 actualyìly: 2019 in the original dataset is a typo

% Detrending using LOESS on an expanding (not rolling) window -> more
% stable, allows for reliable LOESS detrending, better according to O'Brien
% (2023) for performance (NB: for univariate data, it coincides with
% LOWESS)

init_window = 30; % Already enough to get a good number of samples to compute summary statistics, short enough to capture the changepoint (analysed in analysis_Temperature.m)

%% 2 -- Calculate summary statistics

var_matrix_infested = zeros(56,5);
ac_matrix_infested = zeros(56,5);
cv_matrix_infested = zeros(56,5);

var_matrix_uninfested = zeros(56,5);
ac_matrix_uninfested = zeros(56,5);
cv_matrix_uninfested = zeros(56,5);

tt = days(init_window+1: end-1);

% Loop on hives (Infested)
for i = 2 : 6
    % Expanding window
    for j = 1 : (length(days) - init_window - 1)
        y = data_DatasetInternalInfested(1:init_window + j, i);
        y = y(~isnan(y)); % NB: another way to get rid of the NaN (which are anyway just 6) would be to use fillgaps. However, this may induce an autoregressive regularity that is spurious (and knowledgeable of future time points) that would impair the interpretation of the results; it is therefore not used.
        smoothing_infested = smooth(y, 0.5, 'loess'); % 0.5 as in O'Brien and references therein
        detrended_infested = y - smoothing_infested;
        
        % Variance
        var_matrix_infested(j,i-1) = var(detrended_infested);
        % AC1
        coeff = corrcoef(detrended_infested(1:end-1),detrended_infested(2:end));
        ac_matrix_infested(j,i-1) = abs(coeff(1,2));
        %CV
        cv_matrix_infested(j,i-1) = var(detrended_infested) ./ mean(detrended_infested);
    end
end

% Loop on hives (Uninfested)
for i = 2 : 6
    % Expanding window
    for j = 1 : (length(days) - init_window - 1)
        y = data_DatasetInternalUninfested(1:init_window + j, i);
        y = y(~isnan(y));
        smoothing_uninfested = smooth(y, 0.5, 'loess'); % 0.5 as in O'Brien and references therein
        detrended_uninfested = y - smoothing_uninfested;
        
        % Variance
        var_matrix_uninfested(j,i-1) = var(detrended_uninfested);
        % AC1
        coeff = corrcoef(detrended_uninfested(1:end-1),detrended_uninfested(2:end));
        ac_matrix_uninfested(j,i-1) = abs(coeff(1,2));
        %CV
        cv_matrix_uninfested(j,i-1) = var(detrended_uninfested) ./ mean(detrended_uninfested);
    end
end

%% 3 -- Plot indicators

changepoints = datetime(['22-Oct-2018';'18-Sep-2018';'16-Oct-2018';'12-Oct-2018';'07-Oct-2018']); % NB: these are the points where the MEAN temperature already drops - not when the colony collapses. Estimated in analysisTemperature.m
colors = lines(5);

% ------- Infested
% Variance
figure()
hold on
for i = 1 : 5
    if i == 2
        plot(tt(1:21), var_matrix_infested(1:21,i),linewidth=1.5); %earlier stop in collecting data (due to collapse)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    else
        plot(tt, var_matrix_infested(:,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    end
end
ylabel('Var', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='northwest', fontsize=14)

% AC1
figure()
hold on
for i = 1 : 5
    if i == 2
        plot(tt(1:21), ac_matrix_infested(1:21,i),linewidth=1.5); %earlier stop in collecting data (due to collapse)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    else
        plot(tt, ac_matrix_infested(:,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    end
end
ylabel('AC(1)', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='northwest', fontsize=14)


% CV
figure()
hold on
for i = 1 : 5
    if i == 2
        plot(tt(1:21), cv_matrix_infested(1:21,i),linewidth=1.5); %earlier stop in collecting data (due to collapse)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    else
        plot(tt, cv_matrix_infested(:,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
    end
end
ylabel('CV', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='southeast', fontsize=14)



% --------- Uninfested
% Variance
figure()
hold on
for i = 1 : 5
        plot(tt, var_matrix_uninfested(:,i),linewidth=1.5)
end
ylabel('Var', fontsize = 20)
legend(["6P", "7C",  "10Bc",  "VF", "VN"], location='northwest', fontsize=14)

% AC1
figure()
hold on
for i = 1 : 5
        plot(tt, ac_matrix_uninfested(:,i),linewidth=1.5)
end
ylabel('AC(1)', fontsize = 20)
legend(["6P", "7C", "10Bc",  "VF",  "VN"], location='northwest', fontsize=14)

% CV
figure()
hold on
for i = 1 : 5
        plot(tt, cv_matrix_uninfested(:,i),linewidth=1.5)
end
ylabel('CV', fontsize = 20)
legend(["6P",  "7C",  "10Bc",  "VF",  "VN"], location='southwest', fontsize=14)



%% 4 -- ROC analysis
% True positive: a significant signal raised, on the collapsed colonies 
% (infected group minus hive 14 that passed the winter), before the end of
% the sampling or before the colony is decleared dead by colony strength
% count. "Before" = any time after an initial window and up to 1 day before.
% False positive: a significant signal raised on non-collapsed colonies
% Include the mean (it may be just easier to look for that instead of
% relying on CSD EWS...)

tau_p = 0. : 0.01 : 1;  % Threshold to test for ROC

true_pos  = zeros(length(tau_p), 5); % On columns: Var, AC1, CV, Mean, Min
false_pos = zeros(length(tau_p), 5); 

% Var
kendall_matrix_var = zeros(length(var_matrix_uninfested), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(var_matrix_uninfested) -1
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_var(n,m) = corr(var_matrix_infested(1:n, m),var_matrix_infested(2:n+1, m),'type','Kendall');
            if m==2
                if n > 21
                    kendall_matrix_var(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_var(n,m) > tau_p(p))
                     true_pos(p, 1) = true_pos(p,1) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_var(n,m) = corr(var_matrix_uninfested(1:n,m),var_matrix_uninfested(2:n+1,m),'type','Kendall');
            else
                kendall_matrix_var(n,m) = corr(var_matrix_uninfested(1:n,m-5),var_matrix_uninfested(2:n+1,m-5),'type','Kendall');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_var(n,m) > tau_p(p))
                     false_pos(p, 1) = false_pos(p,1) + 1;                 
                end
            end
        end
    end
end

% AC1
kendall_matrix_ac = zeros(length(ac_matrix_uninfested), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(ac_matrix_uninfested) -1
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_ac(n,m) = corr(ac_matrix_infested(1:n, m),ac_matrix_infested(2:n+1, m),'type','Kendall');
            if m==2
                if n > 21
                    kendall_matrix_ac(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_ac(n,m) > tau_p(p))
                     true_pos(p, 2) = true_pos(p,2) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_ac(n,m) = corr(ac_matrix_uninfested(1:n,m),ac_matrix_uninfested(2:n+1,m),'type','Kendall');
            else
                kendall_matrix_ac(n,m) = corr(ac_matrix_uninfested(1:n,m-5),ac_matrix_uninfested(2:n+1,m-5),'type','Kendall');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_ac(n,m) > tau_p(p))
                     false_pos(p, 2) = false_pos(p,2) + 1;                 
                end
            end
        end
    end
end



% CV
kendall_matrix_cv = zeros(length(cv_matrix_uninfested), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(cv_matrix_uninfested) -1
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_cv(n,m) = corr(cv_matrix_infested(1:n, m),cv_matrix_infested(2:n+1, m),'type','Kendall');
            if m==2
                if n > 21
                    kendall_matrix_cv(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_cv(n,m) > tau_p(p))
                     true_pos(p, 3) = true_pos(p,3) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_cv(n,m) = corr(cv_matrix_uninfested(1:n,m),cv_matrix_uninfested(2:n+1,m),'type','Kendall');
            else
                kendall_matrix_cv(n,m) = corr(cv_matrix_uninfested(1:n,m-5),cv_matrix_uninfested(2:n+1,m-5),'type','Kendall');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_cv(n,m) > tau_p(p))
                     false_pos(p, 3) = false_pos(p,3) + 1;                 
                end
            end
        end
    end
end


%Mean
kendall_matrix_mean = zeros(length(ac_matrix_uninfested), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(ac_matrix_uninfested) -1
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_mean(n,m) = corr(data_DatasetInternalInfested(init_window + 1: init_window + n, m),data_DatasetInternalInfested(init_window + 2: init_window + n+1, m),'type','Kendall','rows','complete');
            if m==2
                if n > 21
                    kendall_matrix_mean(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_mean(n,m) > tau_p(p))
                     true_pos(p, 4) = true_pos(p,4) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_mean(n,m) = corr(data_DatasetInternalInfested(init_window + 1: init_window + n, m),data_DatasetInternalInfested(init_window + 2: init_window + n+1, m),'type','Kendall','rows','complete');
            else
                kendall_matrix_mean(n,m) = corr(data_DatasetInternalUninfested(init_window + 1: init_window + n, m-5),data_DatasetInternalUninfested(init_window + 2: init_window + n+1, m-5),'type','Kendall','rows','complete');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_mean(n,m) > tau_p(p))
                     false_pos(p, 4) = false_pos(p,4) + 1;                 
                end
            end
        end
    end
end



%% 5 -- Rates
true_pos_rate = true_pos ./ max(true_pos);  
false_pos_rate = false_pos ./ max(false_pos); 

% AUC
AUC = abs( [trapz(false_pos_rate(1:end, 1), true_pos_rate(1:end, 1)), trapz(false_pos_rate(1:end, 2), true_pos_rate(1:end, 2)), trapz(false_pos_rate(1:end, 3), true_pos_rate(1:end, 3)), trapz(false_pos_rate(1:end, 4), true_pos_rate(1:end, 4))]);

best_tau_p_mean = tau_p( find(round(false_pos_rate(1:end,4),5) == 0.40064, 1) );

x = linspace(0,1);

figure()
hold on
plot(false_pos_rate(1:end, 1),true_pos_rate(1:end, 1),'linewidth',1.5)
plot(false_pos_rate(1:end, 2),true_pos_rate(1:end, 2),'linewidth',1.5)
%plot(false_pos_rate(1:end, 3),true_pos_rate(1:end, 3),'linewidth',1.5) %
%CV is quite useless: it already oscillates like crazy by eye
plot(false_pos_rate(1:end, 4),true_pos_rate(1:end, 4),'linewidth',1.5) 
plot(x,x,'--','color','black')
legend({'Variance','AC(1)','Mean'},'Location','southeast','fontsize',14)
ylabel("Sensitivity",'fontsize',20)
xlabel("1-Specificity",'fontsize',20)
axis square;
hold off
