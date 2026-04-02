%% Extract statistical indicators that may be used as EWS
% Here, for 20-mins sampled temperatures. 

% NB2: when the original data are collated, there is always a sudden drop
% in temperature, likely due to the reset. To avoid spurious results, the
% reset is substituted with an interpolation of adjacent points

clc; clear; close all

%% 1 -- Daily-sampled average temperatures
% Load data
data = readmatrix('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\sonde_2F\all_temperature.xlsx');

data_table = readtable('C:\Users\danie\Documents\POSTDOC_Trento\Projects\Api_bistability\data\sonde_2F\all_temperature.xlsx'); 
days = datetime(data_table.Datetime,'InputFormat', 'dd/MM/yy HH.mm.ss');

% Get min over each day (as in Colin et al., 2025)
min_data = zeros(round((length(data)/72)-1),10);
min_days = days(1:72:end-71,1);
for j=2:11
    for i = 0 : round((length(data)/72))-2
        min_data(i+1,j-1) = min(data(i*71 + 1 : 71 + i*72 , j));
    end
end


% First, I remove the daily trend, then I'll apply LOESS as for the low
% frequency dataset

D12 = LagOp({1 -1},'Lags',[0,71]); % 71 datapoints per day
for i=2:11
    ddata(:,i-1) = filter(D12,data(:,i));
end
len_diff = length(data) - length(ddata);

init_window = 700; % Already enough to get a good number of samples to compute summary statistics, short enough to capture the changepoint (analysed in analysis_Temperature.m)



%% 2 -- Calculate summary statistics

dimension = round((length(ddata) - init_window - 36)/36);

var_matrix_infested = zeros(dimension,5);
ac_matrix_infested = zeros(dimension,5);
cv_matrix_infested = zeros(dimension,5);
spec_matrix_infested = zeros(dimension,5);

var_matrix_uninfested = zeros(dimension,5);
ac_matrix_uninfested = zeros(dimension,5);
cv_matrix_uninfested = zeros(dimension,5);
spec_matrix_uninfested = zeros(dimension,5);

tt = days(init_window +36 : 36 : length(ddata));

% Loop on hives (Infested)
for i = 1 : 5
    % Expanding window
    for j = 1 : (length(tt)) % I do the expanding window twice a day (it's useless to do it every 20 minutes)
        y = ddata(1:init_window + j*36, i);
        y = y(~isnan(y));
        smoothing_infested = smooth(y, 0.5, 'loess'); % 0.5 as in O'Brien and references therein
        detrended_infested = y - smoothing_infested;
        
        % Variance
        var_matrix_infested(j,i) = var(detrended_infested);
        % AC1
        coeff = corrcoef(detrended_infested(1:end-1),detrended_infested(2:end));
        ac_matrix_infested(j,i) = abs(coeff(1,2));
        % CV
        cv_matrix_infested(j,i) = var(detrended_infested) ./ mean(detrended_infested);
        % Spectral reddening
        result = estimate_noise_color(detrended_infested, 0.014);
        spec_matrix_infested(j,i) = result.beta;
    end
end

% Loop on hives (Uninfested)
for i = 6 : 10
    % Expanding window
    for j = 1 : (length(tt))
        y = ddata(1:init_window + j*36, i);
        y = y(~isnan(y));
        smoothing_uninfested = smooth(y, 0.5, 'loess'); % 0.5 as in O'Brien and references therein
        detrended_uninfested = y - smoothing_uninfested;
        
        % Variance
        var_matrix_uninfested(j,i-5) = var(detrended_uninfested);
        % AC1
        coeff = corrcoef(detrended_uninfested(1:end-1),detrended_uninfested(2:end));
        ac_matrix_uninfested(j,i-5) = abs(coeff(1,2));
        % CV
        cv_matrix_uninfested(j,i-5) = var(detrended_uninfested) ./ mean(detrended_uninfested);
        % Spectral reddening
        result = estimate_noise_color(detrended_uninfested, 0.014);
        spec_matrix_uninfested(j,i-5) = result.beta;
    end
end

%% 3 -- Plot indicators
end_vis = 156; % Stop on 5th November as for the other dataset
changepoints = datetime(['22-Oct-2018';'18-Sep-2018';'16-Oct-2018';'12-Oct-2018';'07-Oct-2018']); % NB: these are the points where the MEAN temperature already drops - not when the colony collapses. Estimated in analysisTemperature.m
colors = lines(5);

% ------- Infested
% Variance
figure()
hold on
for i = 1 : 5
    plot(tt(1:end_vis), var_matrix_infested(1:end_vis,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
end
ylabel('Var', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='northwest', fontsize=14)

% AC1
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), ac_matrix_infested(1:end_vis,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
end
ylabel('AC(1)', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='northwest', fontsize=14)


% CV
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), cv_matrix_infested(1:end_vis,i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
end
ylabel('CV', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='southeast', fontsize=14)

% Min
figure()
hold on
for i = 1 : 5
        plot(min_days(1:89), min_data((1:89),i),linewidth=1.5)
        xline(changepoints(i), '--', 'Color',colors(i,:), linewidth=1.5)
end
ylabel('Min', fontsize = 20)
legend(["5N", "", "8R", "", "13R", "", "14R", "", "CR/12", ""], location='southwest', fontsize=14)

% Spectral reddening
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), spec_matrix_infested(1:end_vis,i),linewidth=1.5)
end
ylabel('Spectral slope \beta', fontsize = 20)
legend(["6P", "7C",  "10Bc",  "VF", "VN"], location='southwest', fontsize=14)



% --------- Uninfested
% Variance
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), var_matrix_uninfested(1:end_vis,i),linewidth=1.5)
end
ylabel('Var', fontsize = 20)
legend(["6P", "7C",  "10Bc",  "VF", "VN"], location='northwest', fontsize=14)

% AC1
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), ac_matrix_uninfested(1:end_vis,i),linewidth=1.5)
end
ylabel('AC(1)', fontsize = 20)
legend(["6P", "7C", "10Bc",  "VF",  "VN"], location='northwest', fontsize=14)

% CV
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), cv_matrix_uninfested(1:end_vis,i),linewidth=1.5)
end
ylabel('CV', fontsize = 20)
legend(["6P",  "7C",  "10Bc",  "VF",  "VN"], location='southeast', fontsize=14)

% Min
figure()
hold on
for i = 1 : 5
        plot( min_days(1:89), min_data(1:89,i),linewidth=1.5)
end
ylabel('Min', fontsize = 20)
legend(["6P", "7C",  "10Bc",  "VF", "VN"], location='southwest', fontsize=14)

% Spectral reddening
figure()
hold on
for i = 1 : 5
        plot(tt(1:end_vis), spec_matrix_uninfested(1:end_vis,i),linewidth=1.5)
end
ylabel('Spectral slope \beta', fontsize = 20)
legend(["6P", "7C",  "10Bc",  "VF", "VN"], location='southwest', fontsize=14)


%% 4 -- ROC analysis
% True positive: a significant signal raised, on the collapsed colonies 
% (infected group minus hive 14 that passed the winter), before the end of
% the sampling or before the colony is decleared dead by colony strength
% count. "Before" = any time after an initial window and up to 1 day before.
% False positive: a significant signal raised on non-collapsed colonies
% Include the mean (it may be just easier to look for that instead of
% relying on CSD EWS...)

tau_p = 0. : 0.01 : 1;  % Threshold to test for ROC

true_pos  = zeros(length(tau_p), 6); % On columns: Var, AC1, CV, Mean, Min, Spec Slope
false_pos = zeros(length(tau_p), 6); 

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
            kendall_matrix_mean(n,m) = corr(ddata(init_window + 1: init_window + n, m),ddata(init_window + 2: init_window + n+1, m),'type','Kendall','rows','complete');
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
                kendall_matrix_mean(n,m) = corr(ddata(init_window + 1: init_window + n, m),ddata(init_window + 2: init_window + n+1, m),'type','Kendall','rows','complete');
            else
                kendall_matrix_mean(n,m) = corr(ddata(init_window + 1: init_window + n, m-5),ddata(init_window + 2: init_window + n+1, m-5),'type','Kendall','rows','complete');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_mean(n,m) > tau_p(p))
                     false_pos(p, 4) = false_pos(p,4) + 1;                 
                end
            end
        end
    end
end


%Min
kendall_matrix_min = zeros(length(min_data), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(min_data) - 1 - 10
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_min(n,m) = corr(min_data(10 + 1: 10 + n, m),min_data(10 + 2: 10 + n + 1, m),'type','Kendall','rows','complete');
            if m==2
                if n > 21
                    kendall_matrix_min(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_min(n,m) > tau_p(p))
                     true_pos(p, 5) = true_pos(p,5) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_min(n,m) = corr(min_data(10 + 1: 10 + n, m),min_data(10 + 2: 10 + n+1, m),'type','Kendall','rows','complete');
            else
                kendall_matrix_min(n,m) = corr(min_data(10 + 1: 10 + n, m-5),min_data(10+ 2: 10 + n+1, m-5),'type','Kendall','rows','complete');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_min(n,m) > tau_p(p))
                     false_pos(p, 5) = false_pos(p,5) + 1;                 
                end
            end
        end
    end
end

% Spectral slope
kendall_matrix_spec = zeros(length(spec_matrix_uninfested), 10);  % for each hive, for each point before the end
for m= 1 : 10 % loop on hives
    for n= 1 : length(spec_matrix_uninfested) -1
        if m < 4 || m == 5 % Collapsed colonies (recall that m=4 of the uninfested had not collapsed
            kendall_matrix_spec(n,m) = corr(spec_matrix_infested(1:n, m),spec_matrix_infested(2:n+1, m),'type','Kendall');
            if m==2
                if n > 21
                    kendall_matrix_spec(n,m) = 0;
                end
            end
            for p=1:length(tau_p)
                if (kendall_matrix_spec(n,m) > tau_p(p))
                     true_pos(p, 6) = true_pos(p,6) + 1;                 
                end
            end
        else % Non-collapsed colonies
            if m == 4
                kendall_matrix_spec(n,m) = corr(spec_matrix_uninfested(1:n,m),spec_matrix_uninfested(2:n+1,m),'type','Kendall');
            else
                kendall_matrix_spec(n,m) = corr(spec_matrix_uninfested(1:n,m-5),spec_matrix_uninfested(2:n+1,m-5),'type','Kendall');
            end
            for p=1:length(tau_p)
                if (kendall_matrix_spec(n,m) > tau_p(p))
                     false_pos(p, 6) = false_pos(p,6) + 1;                 
                end
            end
        end
    end
end

%% 5 -- Rates
true_pos_rate = true_pos ./ max(true_pos);  
false_pos_rate = false_pos ./ max(false_pos); 

% AUC
AUC = abs( [trapz(false_pos_rate(1:end, 1), true_pos_rate(1:end, 1)), trapz(false_pos_rate(1:end, 2), true_pos_rate(1:end, 2)), trapz(false_pos_rate(1:end, 3), true_pos_rate(1:end, 3)), trapz(false_pos_rate(1:end, 4), true_pos_rate(1:end, 4)), trapz(false_pos_rate(1:end, 5), true_pos_rate(1:end, 5)), trapz(false_pos_rate(1:end, 6), true_pos_rate(1:end, 6)) ]);

best_tau_p_mean = tau_p( find(round(false_pos_rate(1:end,4),5) == 0.36402, 1) );

x = linspace(0,1);

figure()
hold on
plot(false_pos_rate(1:end, 1),true_pos_rate(1:end, 1),'linewidth',1.5)
plot(false_pos_rate(1:end, 2),true_pos_rate(1:end, 2),'linewidth',1.5)
%plot(false_pos_rate(1:end, 3),true_pos_rate(1:end, 3),'linewidth',1.5) %
%CV is quite useless: it already oscillates like crazy by eye
plot(false_pos_rate(1:end, 4),true_pos_rate(1:end, 4),'linewidth',1.5) 
plot(false_pos_rate(1:end, 5),true_pos_rate(1:end, 5),'linewidth',1.5) 
plot(false_pos_rate(1:end, 6),true_pos_rate(1:end, 6),'linewidth',1.5) 
plot(x,x,'--','color','black')
legend({'Variance','AC(1)','Mean','Min', '\beta'},'Location','southeast','fontsize',15)
ylabel("Sensitivity",'fontsize',20)
xlabel("1-Specificity",'fontsize',20)
axis square;
hold off




%% 6 -- Information about noise (before and after the changepoint)
chpts = [5294, 2828, 4844, 4556, 4196, 4845, 4700, 4556, 4340, 4916]; %indeces of changepoints


for i = 1:length(chpts)
    SG_filtered_before = sgolayfilt(data(1:chpts(i), i+1),2,35); % data, order=2, window size = 35 (half a day)
    z = data(1:chpts(i), i+1) - SG_filtered_before;
    SG_detrended_before{i} = z(~isnan(z));
        
    SG_filtered_after = sgolayfilt(data(chpts(i):end, i+1),2,35); % data, order=2, window size = 35 (half a day)
    zz = data(chpts(i):end, i+1) - SG_filtered_after;
    SG_detrended_after{i} = zz(~isnan(zz));
end

% --- Noise properties

% % Before the changepoint
% Distribution: close to Gaussian?
apiaries = ["5N",  "8R",  "13R",  "14R",  "CR/12", "6P", "7C",  "10Bc",  "VF", "VN"];
for i =1 :10
    figure()
    histogram(SG_detrended_before{i}, 'Normalization', 'pdf') 
    txt = ['\mu=  ' num2str(skewness(SG_detrended_before{i}, 0))];
    text(1,2.5,txt, fontsize=15)
    txt = ['\gamma^*=  ' num2str(kurtosis(SG_detrended_before{i}, 0)-3)]; % excess of kurtosis
    text(1,2,txt, fontsize=15)
    title(apiaries(i), fontsize=12)
end

% Frequency properties: which color?
for i =1 :10
    result = estimate_noise_color(SG_detrended_before{i}, 0.014);
    fprintf('Beta = %.2f → %s\n', result.beta, result.color);
    plot_psd_fit(result, apiaries(i));
end

% % After the changepoint
% Distribution: close to Gaussian?
apiaries = ["5N",  "8R",  "13R",  "14R",  "CR/12", "6P", "7C",  "10Bc",  "VF", "VN"];
for i =1 :10
    figure()
    histogram(SG_detrended_after{i}, 'Normalization', 'pdf') 
    txt = ['\mu =  ' num2str(skewness(SG_detrended_after{i}, 0))];
    text(1,1,txt, fontsize=15)
    txt = ['\gamma^* =  ' num2str(kurtosis(SG_detrended_after{i}, 0)-3)]; % excess of kurtosis
    text(1,0.5,txt, fontsize=15)
    %txt = ['Var =  ' num2str(var(SG_detrended_after{i}, 0))]; % excess of kurtosis
    %text(1,0.1,txt, fontsize=15)
    title(apiaries(i), fontsize=12)
end

% --- Frequency properties: which color?
for i =1 :10
    result = estimate_noise_color(SG_detrended_after{i}, 0.014);
    fprintf('Beta = %.2f → %s\n', result.beta, result.color);
    plot_psd_fit(result, apiaries(i));
end