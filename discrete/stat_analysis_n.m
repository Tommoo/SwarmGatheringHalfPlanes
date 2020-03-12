clear; close all; clc;

square_size = 50;

nb_trials = 1000; %Number of runs of the algorithms (nb of rand initialisations and then run the algo on it)
dt = 1;%Time steps of resampling the angles
v = 1; %Speed: Step size of v*dt in the discrete case

plot_hist_evolution = false;
plot_mean_evolution = false;

save_res_fig = true;



DIR_RES = 'res/';
DIR_RES_CONV_TIME = [DIR_RES,'conv_time/'];

files = dir([DIR_RES_CONV_TIME,'conv_time_*_square_size_',num2str(square_size),'_*.txt']);
files_names = {files.name};

t_all = zeros(nb_trials,length(files_names));
t_mean = zeros(1,length(files_names));
n = zeros(1,length(files_names));
std = zeros(1,length(files_names));
t_cumsum = zeros(nb_trials,length(files_names));
t_distrib = cell(1,length(files_names));

%Sort with increasing n to avoid bad surprises
for file_idx = 1:length(files_names)
    fname = files_names{file_idx};
    
    str = strsplit(fname,'_');
    n(file_idx) = str2num(str{4});
end
[n,sort_idx] = sort(n,'ascend'); 


legends_mean_cv_all_normalised = {}; %Normalise in same plot cumsum cv time, by the last mean cv time, so that we can see in a same plot
h_mean_cv_all_normalised = figure('visible','off'); hold on

list_h_hist_all_normalised = cell(1,length(files_names)); %Lists the hFig pointer to the figures for each file of the histogram summary


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
            title(['n = ',num2str(n(file_idx)),', nb_{trials} = ',num2str(i)])
            xlabel('Time of convergence')
            drawnow
            %hold off
        end
    end
    
    if plot_mean_evolution
        figure;
        plot(1:nb_trials,t_cumsum(:,file_idx)./(1:nb_trials)');
        title([convertCharsToStrings(['Average time of convergence Sum_{i=1..k}({(t_{cv})}_i/k)   for k\in [1,',num2str(nb_trials),']']),...
               convertCharsToStrings(['with n = ',num2str(n(file_idx))])])
        xlabel('Number of trials');
        ylabel('Average time of convergence')
    end
    
    %Mean evolution all in one figure
    %figure(h_mean_cv_all_normalised); %Visible at each step
    set(0,'CurrentFigure',h_mean_cv_all_normalised); %Invisible: calls the current figure to be h_mean_cv_all_normalised
    plot(1:nb_trials,(t_cumsum(:,file_idx)./(1:nb_trials)')/mean(t_conv),'LineWidth',1.5);
    legends_mean_cv_all_normalised{end+1} = ['n = ',num2str(n(file_idx))];
    
    %Histogram evolution all in one figure
    legends_hist_all_normalised = {}; %Normalise such that sum all mass equals to 1 for plotting on same figure (since we add more an more trials)
    h_hist_all_normalised = figure('visible','off');
    hFig = histogram(t_conv,'Normalization','probability'); %Preprocessing to get the hist edges shared by everyone
    hist_edges = hFig.BinEdges;
    clf; hold on %Remove the first histogram computation
    for i = [50,100,300,500,1000] %Only 10 histograms
        t_conv_sub = t_conv(1:i);
        hFig = histogram(t_conv_sub,hist_edges,'Normalization','probability','facealpha',0.2);
        legends_hist_all_normalised{end+1} = [num2str(i),' trials'];
    end
    hFig.EdgeColor = 'r'; %For the last histogram
    hFig.LineWidth = 1.5; %For the last histogram
    title(['Normalised distribution of convergence time with n = ',num2str(n(file_idx)),' agents'])
    xlabel('Time of convergence')
    ylabel(["Normalised count","Count_{bin}/Count_{total}"])
    legend(legends_hist_all_normalised)
    list_h_hist_all_normalised{file_idx} = h_hist_all_normalised;
end

set(h_mean_cv_all_normalised,'visible','on')
figure(h_mean_cv_all_normalised);
legend(legends_mean_cv_all_normalised);
xlabel('Number of trials','FontWeight','bold');
%ylabel(["Normalised mean convergence time",convertCharsToStrings(['Sum_{i=1..k}({(t_{cv})}_i/k)/Sum_{i=1..1000}({(t_{cv})}_i/1000)  for  k\in [1,',num2str(nb_trials),']'])],'FontWeight','bold')
ylabel(["Normalised mean convergence time",convertCharsToStrings(['Mean_{k}(t_{cv})/Mean_{1000}(t_{cv})  for  k\in [1,',num2str(nb_trials),']'])],'FontWeight','bold')
title('Normalised mean convergence time')

if save_res_fig 
    saveas(gcf,['res/conv_time/figures/discrete_nAnalysis_cumsum_tcv_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/discrete_nAnalysis_cumsum_tcv_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
end

for file_idx = 1:length(files_names)
    set(list_h_hist_all_normalised{file_idx},'visible','on');
    
    if save_res_fig
        saveas(gcf,['res/conv_time/figures/discrete_nAnalysis_hist_tcv_n_',num2str(n(file_idx)),'_square_size_',num2str(square_size),'_v_',num2str(v)],'fig')
        saveas(gcf,['res/conv_time/figures/discrete_nAnalysis_hist_tcv_n_',num2str(n(file_idx)),'_square_size_',num2str(square_size),'_v_',num2str(v)],'epsc')
        %Eps save does not handle transparency of histogram so it rasterises instead :-( --> Go for pdf
        
        %Save figure in A4 page (bad: lots of white space)
        %saveas(gcf,['res/conv_time/figures/discrete_nAnalysis_hist_tcv_n_',num2str(n(file_idx)),'_square_size_',num2str(square_size),'_v_',num2str(v)],'pdf')
        %Save figure with minimum white space
        %https://fr.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,['res/conv_time/figures/discrete_nAnalysis_hist_tcv_n_',num2str(n(file_idx)),'_square_size_',num2str(square_size),'_v_',num2str(v)],'-dpdf')
    end
end



figure; hold on
plot(n,t_mean,'+-')
xlabel('n')
ylabel('Mean convergence time')
legend(['square size = ',num2str(square_size)])

figure;
errorbar(n,t_mean,std)
xlabel('n')
ylabel('Mean convergence time')
legend(['square size = ',num2str(square_size)])


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
legend('Results','Empirical mean and std','Weighted LS regression')
%title(['Convergence time (Initial square size = ',num2str(square_size),'\times',num2str(square_size),')'])
title('Convergence time')

if save_res_fig
    saveas(gcf,['res/conv_time/figures/discrete_square_size_',num2str(square_size),'_dt_',num2str(dt),'_v_',num2str(v)],'fig')
    saveas(gcf,['res/conv_time/figures/discrete_square_size_',num2str(square_size),'_dt_',num2str(dt),'_v_',num2str(v)],'epsc')
end

