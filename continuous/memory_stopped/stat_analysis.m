clear; close all; clc;


n = 10; %nb agents
nb_trials = 1000; %Number of runs of the algorithms (nb of rand initialisations and then run the algo on it)

plot_hist_evolution = false;
plot_mean_evolution = false;

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

files = dir([DIR_RES_CONV_TIME,'conv_time_n_',num2str(n),'*.txt']);
files_names = {files.name};

t_all = zeros(nb_trials,length(files_names));
t_mean = zeros(1,length(files_names));
delta = zeros(1,length(files_names));
std = zeros(1,length(files_names));
t_cumsum = zeros(nb_trials,length(files_names));
t_distrib = cell(1,length(files_names));

%Sort with increasing delta to avoid bad surprises
for file_idx = 1:length(files_names)
    fname = files_names{file_idx};
    
    str = strsplit(fname,'_');
    delta(file_idx) = str2num(str{6});
end
[delta,sort_idx] = sort(delta,'ascend'); 


for file_idx = 1:length(files_names)
    fname = files_names{sort_idx(file_idx)};
    
    %str = strsplit(fname,'_');
    %delta(file_idx) = str2num(str{6});
    
    
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
            title(['n = ',num2str(n),', delta = ',num2str(delta(file_idx)),', nb_{trials} = ',num2str(i)])
            xlabel('Time of convergence')
            drawnow
            %hold off
        end
    end
    
    if plot_mean_evolution
        figure;
        plot(1:nb_trials,t_cumsum(:,file_idx)./(1:nb_trials)');
        title(['n = ',num2str(n),', delta = ',num2str(delta(file_idx))])
        xlabel(['Sum_{i=1..k}(T_i/k) for k=1..nb_{trials}=1..',num2str(nb_trials)]);
        ylabel('Average time of convergence')
    end
end


bound = 4*n*ceil(n*(n-1)*sqrt(2)./(2*delta*(1-sqrt(1+(tan(pi/(4*n)))^2-2*(tan(pi/(4*n)))*sin(pi/(2*n))))));

figure; hold on
plot(delta,t_mean,'+-')
xlabel('delta')
ylabel('Mean time of convergence')
legend(['n = ',num2str(n)])

figure; hold on
plot(delta,bound,'r+-');
plot(delta,t_mean,'+-')
xlabel('delta')
ylabel('Mean time of convergence')
legend(['Worst case bound'],['n = ',num2str(n)])

figure;
plot(1./delta,t_mean,'+-')
xlabel('1/delta')
ylabel('Mean time of convergence')
legend(['n = ',num2str(n)])

figure; hold on
plot(1./delta,bound,'r+-')
plot(1./delta,t_mean,'+-')
xlabel('1/delta')
ylabel('Mean time of convergence')
legend('Worst case bound',['n = ',num2str(n)])

figure;
errorbar(1./delta,t_mean,std)
xlabel('1/delta')
ylabel('Mean time of convergence')
legend(['n = ',num2str(n)])


