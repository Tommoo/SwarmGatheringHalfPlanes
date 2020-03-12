clear; close all; clc;

warning('off','all')
warning

n = 10; %nb agents
nb_steps = realmax;%1000; %nb of steps
dt = 1;%0.05; %Time steps of resampling the angles
v = 1; %Speed: Step size of v*dt in the discrete case
ddt = dt/200;%dt/100; %Discretisation of the continuous time %We do not resample alpha but do test half plane
              %Mini step of v*ddt
square_size = 1; %Size of the square for rand uniform initalisation [0,square_size]x[0,square_size]
nb_trials = 1000; %Number of runs of the algorithms (nb of rand initialisations and then run the algo on it)

%delta_blind = nan; %Agents are not blind to close agents
delta_blind = 0.00625;%0.01; %If the agents are too close they don't see each other (otherwise getting too close
                    %will force the agents to not be able to move any more)

limited_range = nan; %limited range. Nan if none
%%limited_range = 0.15; %limited range. Nan if none

compute_circle = true; %Computes the enclosing circle
early_stop_convergence = compute_circle & true; %Stop if we have converged in a circle of radius delta_blind

disp_realTime = false; %Shows real time progress
save_video = false; %Saves the video

save_conv_time_to_file = compute_circle & true; %Save the first time to convergence in txt file

FPS = 10; %Nb frames per second for saved video

%rng(42); %Random seed


if save_conv_time_to_file
    fname_conv = ['res/conv_time/conv_time_n_',num2str(n),'_delta_',num2str(delta_blind),'_square_size_',num2str(square_size),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'_ddt_',num2str(ddt),'.txt'];
    fileID = fopen(fname_conv,'a');
end


%%%%%%%%% %For parfor
t_conv = zeros(nb_trials,1);
%%%%%%%%%


parfor trial = 1:nb_trials
    fprintf("Trial %d/%d...",trial,nb_trials);
    
    tic
    
    if save_video
        if isnan(limited_range)
            outputVideo = VideoWriter(['ants_randomHalfPlaneSensor_memoryMoved_continuous_not_limited_n_',num2str(n),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'_ddt_',num2str(ddt),'_delta_blind_',num2str(delta_blind),'_trial_',num2str(nb_trials),'.avi']);
        else
            outputVideo = VideoWriter(['ants_randomHalfPlaneSensor_memoryMoved_continuous_limited_range_',num2str(limited_range),'_n_',num2str(n),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'_ddt_',num2str(ddt),'_delta_blind_',num2str(delta_blind),'_trial_',num2str(nb_trials),'.avi']);
        end
        outputVideo.FrameRate = FPS;
        open(outputVideo)
    end
    
    
    
    
    p = rand(2,n)*square_size; %Initial positions
    
    scatter(p(1,:),p(2,:));
    axis equal
    xlim([0,square_size])
    ylim([0,square_size])
    
    
    
    for step = 1:nb_steps
        %fprintf("step = %d\n",step);
        alpha = rand(1,n)*2*pi; %Random direction: the opposite of the movement
        stopped = false(1,n); %Boolean memory of each agent whether it stoped in the step
        
        
        for step_mini = 1:dt/ddt
            p_next = p;
            for i = 1:n
                dir_i = p(:,i) - [p(:,1:i-1),p(:,i+1:end)]; %Vectors to other agents
                dists_i = sqrt(sum(dir_i.^2,1)); %Distances to other agents
                u_i = dir_i./dists_i; %Unit vectors to other agents
                
                u_alpha = [cos(alpha(i));sin(alpha(i))]; %Opposite of potential speed vector of agent i
                v_i = -u_alpha*v;
                
                scal_i = sum(u_alpha.*u_i,1); %Scalar products of u_alpha and u_i
                
                visible_i = [];
                if isnan(limited_range) %No limit to range
                    visible_i = scal_i <= 0;
                elseif limited_range>0 %Finite limit to range
                    visible_i = scal_i <= 0 & dists_i <= limited_range;
                else
                    error("The parameter range should be either nan or strictly positive");
                end
                
                if isnan(delta_blind) %No blind zone
                    visible_i = visible_i; %Don't change anything
                elseif delta_blind>0 %Close blind zone
                    visible_i = visible_i & dists_i >= delta_blind;
                else
                    
                end
                
                if sum(visible_i) > 0 %If there is one point that is visible then we don't do anything
                    v_i = [0;0];
                    stopped(i) = true;
                elseif stopped(i) == true %If the point had previously stoped then it cannot move again
                    v_i = [0;0];
                end
                
                p_next(:,i) = p(:,i) + v_i*ddt;
            end
            
            p = p_next;
            
            
            
            %{
            if disp_realTime
                scatter(p(1,:),p(2,:));
                axis equal
                xlim([0,square_size])
                ylim([0,square_size])
                drawnow
            end
            %}
            
        end
        
        if compute_circle
            [x_c,y_c,r] = SmallestEnclosingCircle(p(1,:),p(2,:));
        end
        %{
        if disp_realTime
            scatter(p(1,:),p(2,:));
            axis equal
            
            if compute_circle
                hold on
                viscircles([x_c,y_c],r,'Color','r');
                plot(x_c,y_c,'r+')
                hold off
            end
            
            xlim([0,square_size])
            ylim([0,square_size])
            drawnow
        end
        
        
        
        if save_video
            F = getframe(gcf);
            [X, Map] = frame2im(F);
            writeVideo(outputVideo,X);
        end
        %}
        if early_stop_convergence
            if r<=delta_blind
                if save_conv_time_to_file
                    %fprintf(fileID, '%d\n', step);
                    t_conv(trial) = step;
                end
                break
            end
        end
    end
    
    
    %{
    scatter(p(1,:),p(2,:));
    axis equal
    xlim([0,square_size])
    ylim([0,square_size])
    if compute_circle
        hold on
        viscircles([x_c,y_c],r,'Color','r');
        plot(x_c,y_c,'r+')
        hold off
    end
    %}
    if save_video
        close(outputVideo);
    end
    
    time = toc;
    fprintf(" %.2fs\n",time);
end


%%%%%%%%% %For parfor
if save_conv_time_to_file
    for i=1:nb_trials
        fprintf(fileID, '%d\n', t_conv(i));
    end
end
%%%%%%%%%


if save_conv_time_to_file
   fclose(fileID);
end

