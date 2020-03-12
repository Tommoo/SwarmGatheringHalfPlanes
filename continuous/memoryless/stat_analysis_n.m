clear; close all; clc;


delta = 0.02; %nb agents
nb_trials = 1000; %Number of runs of the algorithms (nb of rand initialisations and then run the algo on it)
square_size = 1;
v = 1; %Speed: Small step size of v*ddt (with ddt<=dt=1)

plot_hist_evolution = false;
plot_mean_evolution = false;

save_res_fig = false;

%{
nb_steps = realmax;%1000; %nb of steps
dt = 1;%0.05; %Time steps of resampling the angles
v = 1; %Speed: Step size of v*dt in the discrete case
ddt = dt/500;%dt/100; %Discretisation of the continuous time %We do not resample alpha but do test half plane
              %Mini step of v*ddt
square_size = 1; %Size of the square for rand uniform initalisation [0,square_size]x[0,square_size]

%delta_blind = nan; %Agents are not blind to close agents
delta_blind = 0.005;%0.01; %If the agents are too close they don't see each other (otherwise getting too close
                    %will force the agents to not be able to move any more)
%}
DIR_RES = 'res/';
DIR_RES_CONV_TIME = [DIR_RES,'conv_time/'];

files = dir([DIR_RES_CONV_TIME,'conv_time_*_delta_',num2str(delta),'_*.txt']);
files_names = {files.name};

t_all = zeros(nb_trials,length(files_names));
t_mean = zeros(1,length(files_names));
n = zeros(1,length(files_names));
std = zeros(1,length(files_names));
t_cumsum = zeros(nb_trials,length(files_names));
t_distrib = cell(1,length(files_names));

%Sort with increasing delta to avoid bad surprises
for file_idx = 1:length(files_names)
    fname = files_names{file_idx};
    
    str = strsplit(fname,'_');
    n(file_idx) = str2num(str{4});
end
[n,sort_idx] = sort(n,'ascend'); 

legends_mean_cv_all_normalised = {}; %Normalise in same plot cumsum cv time, by the last mean cv time, so that we can see in a same plot
h_mean_cv_all_normalised = figure('visible','off'); hold on

for file_idx = 1:length(files_names)
    fname = files_names{sort_idx(file_idx)};
    
    %str = strsplit(fname,'_');
    %n(file_idx) = str2num(str{4});
    
    
    fileID = fopen([DIR_RES_CONV_TIME,fname],'r');
    t_conv = fscanf(fileID,'%d');
    fclose(fileID);
    
    t_all(:,file_idx) = t_conv;

    t_mean(file_idx) = mean(t_conv); %Mean
    std(file_idx) = sqrt((1/(nb_trials-1))*sum((t_conv-mean(t_conv)).^2)); %Standard deviation (unbiased)
    
    t_cumsum(:,file_idx) = cumsum(t_conv);
    
    if plot_hist_evolution
        t_distrib_array = zeros(nb_trials,nb_trials);
        figure;
        for i = 2:nb_trials
            t_conv_sub = t_conv(1:i);
            %t_distrib_array = sqrt(i)*(t_conv_sub - mean(t_conv_sub))/(sqrt((1/(i-1))*sum((t_conv_sub-mean(t_conv_sub)).^2)));
            %t_distrib_array = (t_conv_sub - mean(t_conv_sub))/(sqrt((1/(i-1))*sum((t_conv_sub-mean(t_conv_sub)).^2)));
            t_distrib_array = t_conv_sub;
            
            t_distrib{file_idx} = t_distrib_array;
            
            histogram(t_distrib_array,'Normalization','probability');
            %histfit(t_distrib_array);
            %hold on;
            %x = linspace(min(t_distrib_array),max(t_distrib_array),1000);
            %plot(x,normpdf(x,0,1));
            title(['n = ',num2str(n),', delta = ',num2str(delta),', nb_{trials} = ',num2str(i)])
            xlabel('Time of convergence')
            drawnow
            %hold off
        end
    end
    
    if plot_mean_evolution
        figure;
        plot(1:nb_trials,t_cumsum(:,file_idx)./(1:nb_trials)');
        title([convertCharsToStrings(['Average time of convergence Sum_{i=1..k}({(t_{cv})}_i/k)   for k\in [1,',num2str(nb_trials),']']),...
               convertCharsToStrings(['with n = ',num2str(n),', delta = ',num2str(delta)])])
        xlabel('Number of trials');
        title(['n = ',num2str(n),', delta = ',num2str(delta)])
        ylabel('Average time of convergence')
    end
    
    figure(h_mean_cv_all_normalised);
    plot(1:nb_trials,(t_cumsum(:,file_idx)./(1:nb_trials)')/mean(t_conv),'LineWidth',1.5);
    legends_mean_cv_all_normalised{end+1} = ['n = ',num2str(n(file_idx))];
end

set(h_mean_cv_all_normalised,'visible','on')
figure(h_mean_cv_all_normalised);
legend(legends_mean_cv_all_normalised);
xlabel('Number of trials','FontWeight','bold');
%ylabel(["Normalised mean convergence time",convertCharsToStrings(['Sum_{i=1..k}({(t_{cv})}_i/k)/Sum_{i=1..1000}({(t_{cv})}_i/1000)  for  k\in [1,',num2str(nb_trials),']'])],'FontWeight','bold')
ylabel(["Normalised mean convergence time",convertCharsToStrings(['Mean_{k}(t_{cv})/Mean_{1000}(t_{cv})  for  k\in [1,',num2str(nb_trials),']'])],'FontWeight','bold')
title(['Normalised mean convergence time with \delta = ',num2str(delta)])

if save_res_fig
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_cumsum_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_cumsum_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
end


bound = 8*n.*ceil(n.*(n-1).*sqrt(2)./(2*delta.*(1-sqrt(1+(tan(pi./(4*n))).^2-2*(tan(pi./(4*n))).*sin(pi./(2*n))))));
bound_all_plot = 8*(1:n(end)).*ceil((1:n(end)).*((1:n(end))-1).*sqrt(2)./(2*delta.*(1-sqrt(1+(tan(pi./(4*(1:n(end))))).^2-2*(tan(pi./(4*(1:n(end))))).*sin(pi./(2*(1:n(end)))))))); %For plotting smoothly

figure; hold on
plot(n,t_mean,'+-')
xlabel('n')
ylabel('Mean convergence time')
legend(['delta = ',num2str(delta)])

figure; hold on
plot(n,bound,'r+-');
plot(n,t_mean,'+-')
xlabel('n')
ylabel('Mean convergence time')
legend(['Worst case bound'],['delta = ',num2str(delta)])


figure;
errorbar(n,t_mean,std)
xlabel('n')
ylabel('Mean convergence time')
legend(['delta = ',num2str(delta)])



%Crude res figure (delta,t_cv)
h = figure; hold on        
set(h,'Position',[360,333,560,365])
n_all = repmat(n,1000,1);
scatter(n_all(:),t_all(:),'.');
errorbar(n,t_mean,std,'+-','LineWidth',2,'MarkerSize',3);
xlabel('n','FontWeight','bold')
ylabel('Mean convergence time','FontWeight','bold')
lgd = legend('Results','Empirical mean and std');
set(lgd,'Position',[0.61 0.8247 0.2429 0.0726])
%title(['Empirical convergence time (Initial square size = ',num2str(square_size),'\times',num2str(square_size),')'])
title('Empirical convergence time')

if save_res_fig
    %Need to convert . in 0.02 to - as 0-02 as otherwise matlab thinks it's not a fig file type 
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_crude_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_crude_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
end

%Bound figure (n,t_cv)
h = figure; hold on        
set(h,'Position',[360,333,560,365])
plot(n,t_mean,'+-','LineWidth',1.5,'MarkerSize',3)
l = line(); %Dirty workaround to get empty non plotted line but correct legend with line and marker
l.XData = [];
l.YData = [];
l.Color = 'r';
l.Marker = '+';
l.LineWidth = 1.5;
l.MarkerSize = 3;
plot(n,bound,'r+','LineWidth',1.5,'MarkerSize',3); %Without line (not smooth)
plot(1:n(end),bound_all_plot,'r-','LineWidth',1.5,'MarkerSize',3) %Just the smooth line (no marker)
xlabel('n','FontWeight','bold')
ylabel('Mean convergence time','FontWeight','bold')
lgd = legend('Empirical mean and std','Theoretical bound');
set(lgd,'Position',[0.61 0.8247 0.2429 0.0726])
%title(['Theoretical bound on the convergence time (Initial square size = ',num2str(square_size),'\times',num2str(square_size),')'])
title('Theoretical bound on the convergence time')

if save_res_fig
    %Need to convert . in 0.02 to - as 0-02 as otherwise matlab thinks it's not a fig file type 
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_bound_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_bound_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
end


%Interpolation figure (n,t_cv)
h = figure; hold on        
set(h,'Position',[360,333,560,365])
n_all = repmat(n,1000,1);
scatter(n_all(:),t_all(:),'.');
errorbar(n,t_mean,std,'+','LineWidth',2,'MarkerSize',3);
fitline = fit(n',t_mean','poly1','Weights',1./(std'.^2));
hFig = plot(fitline,'k-');
set(hFig, 'LineWidth',1.5);
xlabel('n','FontWeight','bold')
ylabel('Mean convergence time','FontWeight','bold')
lgd = legend('Results','Empirical mean and std','Weighted LS regression');
set(lgd,'Position',[0.61 0.8247 0.2429 0.0726])
%title(['Convergence time (Initial square size = ',num2str(square_size),'\times',num2str(square_size),')'])
title('Convergence time')
ylim([0,4000])

if save_res_fig
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/continuous_nAnalysis_n_tcv_delta_',strrep(num2str(delta),'.','-'),'_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
end
