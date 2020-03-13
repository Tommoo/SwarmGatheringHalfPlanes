clear; close all; clc;

n = 5; %nb agents
nb_steps = realmax;%4000; %nb of steps
dt = 1;
v = 1; %Speed: Step size of v*dt
square_size = 50; %Size of the square for rand uniform initalisation [0,square_size]x[0,square_size]
nb_trials = 1000; %Number of runs of the algorithms (nb of rand initialisations and then run the algo on it)

limited_range = nan; %limited range. Nan if none
%limited_range = 0.1; %limited range. Nan if none

compute_circle = true; %Computes the enclosing circle
early_stop_convergence = compute_circle & true; %Stop if we have converged in a circle of radius delta_blind

disp_realTime = false; %Shows real time progress
save_video = false; %Saves the video

save_conv_time_to_file = compute_circle & false; %Save the first time to convergence in txt file

%Memory of position or r display is not compatible with parfoor loop
nb_steps_mod_memory = 50; %We memorise the position of the agents on the convex hull every nb_steps_mod_memory steps (starting at 1)
memorise_pos_display = (nb_steps < realmax) & false; %Memorises the position of the agents on the convex hull ever nb_steps_mod_memory steps (for display purposes)
memorise_r_display = (nb_steps < realmax) & compute_circle & false; %Memorises the radius of the enclosing circle of the agents at every step (for display purposes)

FPS = 10; %Nb frames per second for saved video

%seed = 42;
%rng(seed); %Random seed


if save_conv_time_to_file
    fname_conv = ['res/conv_time/conv_time_n_',num2str(n),'_square_size_',num2str(square_size),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'.txt'];
    fileID = fopen(fname_conv,'a');
end


%%%%%%%%% %For parfor
t_conv = zeros(nb_trials,1);
%%%%%%%%%


for trial = 1:nb_trials %Compatible with memory of pos and r display
%parfor trial = 1:nb_trials %Not compatible with memory of pos and r display
    fprintf("Trial %d/%d...",trial,nb_trials);
    
    tic
    
    if save_video
        if isnan(limited_range)
            outputVideo = VideoWriter(['ants_randomHalfPlaneSensor_memoryNone_discrete_not_limited_n_',num2str(n),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'.avi']);
        else
            outputVideo = VideoWriter(['ants_randomHalfPlaneSensor_memoryNone_discrete_limited_range_',num2str(limited_range),'_n_',num2str(n),'_nb_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'.avi']);
        end
        outputVideo.FrameRate = FPS;
        open(outputVideo)
    end
    
    p_0 = rand(2,n)*square_size; %Initial positions
    
    p = p_0;
    
    scatter(p(1,:),p(2,:));
    axis equal
    xlim([0,square_size])
    ylim([0,square_size])
    
    
    if memorise_pos_display
        nb_memories = floor(nb_steps/nb_steps_mod_memory)+1; %+1 to keep the initial one and after the n updates
        p_memories = cell(1,nb_memories+1);  
        for memory_idx = 1:nb_memories
            p_memories{memory_idx} = nan(2,n);
        end
        p_memories{1} = p;
    end
    
    if memorise_r_display
        r_memories_all = nan(1,nb_steps+1); %+1 since we want the initial r and after the the nb_steps updates of r (thus nb_steps+1)
        [~,~,r] = SmallestEnclosingCircle(p(1,:),p(2,:));
        r_memories_all(1) = r;
    end
    
    
    for step = 1:nb_steps
        %fprintf("step = %d\n",step);
        %if mod(step,100) ==0
        %    fprintf("\t step = %d\n",step);
        %end
        alpha = rand(1,n)*2*pi; %Random direction: the opposite of the movement
        
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
            
            if sum(visible_i) > 0 %If there is one point that is visible then we don't do anything
                v_i = [0;0];
            end
            
            p_next(:,i) = p(:,i) + v_i*dt;
        end
        
        p = p_next;
        
        
        if compute_circle
            [x_c,y_c,r] = SmallestEnclosingCircle(p(1,:),p(2,:));
        end
        
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
        
        
        if memorise_pos_display && mod(step-1,nb_steps_mod_memory) == 0 %Time to take a screenshot of convex hull, there has been step-1 updates in this screenshot
            p_memories{1+1+(step-1)/nb_steps_mod_memory} = p; %Extra +1 since the first value is for the initial condition (the second +1 since indexing starts at 1)
        end
        
        if memorise_r_display
            r_memories_all(step+1) = r; %Extra +1 since the first value is for the initialisation
        end
        
        if save_video
            F = getframe(gcf);
            [X, Map] = frame2im(F);
            writeVideo(outputVideo,X);
        end
        
        if early_stop_convergence
            if r<=v*dt
                if save_conv_time_to_file
                    %fprintf(fileID, '%d\n', step); %Not great for parfor
                    t_conv(trial) = step;           %Suitable for parfor
                end
                break
            end
        end
    end
    

    
    
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
    
    if save_video
        close(outputVideo);
    end
    
    time = toc;
    fprintf(" %.2fs\n",time);
end


if memorise_pos_display
    hFig = figure; 
    if memorise_r_display
        set(hFig,'Position',[100,100,770,420])
        subplot(1,3,2:3)
    end
    hold on
    plot(p_0(1,:),p_0(2,:),'kx');
    for memory_idx = 1:nb_memories
        p_memory = p_memories{memory_idx};
        if sum(isnan(p_memory))>0
            break %We stopped the loop before the max number of iterations, thus we do not have a screenshot, no need to continue
        else
            conv_hull_ind = convhull(p_memory(1,:),p_memory(2,:));
            plot(p_memories{memory_idx}(1,conv_hull_ind),p_memories{memory_idx}(2,conv_hull_ind),'b:o')
        end
    end
    axis equal    
    xlim([0,square_size])
    ylim([0,square_size])
    title('Ants'' dynamics')
end

if memorise_r_display
    if memorise_pos_display
        subplot(1,3,1)
        hold on
    else
        figure; hold on
    end
    if sum(isnan(r_memories_all)) > 0
        r_mem_idx_nan = find(isnan(r_memories_all),1); %First index of nan
        plot(0:r_mem_idx_nan-1-1,r_memories_all(1:r_mem_idx_nan-1),'-o','MarkerSize',1);
        plot(0:r_mem_idx_nan-1-1,ones(1,r_mem_idx_nan-1),'r--')
        plot(0:r_mem_idx_nan-1-1,2*ones(1,r_mem_idx_nan-1),'k-')
    else
        plot(0:nb_steps,r_memories_all,'-o','MarkerSize',1)
        plot(0:nb_steps,ones(1,nb_steps+1),'r--')
        plot(0:nb_steps,2*ones(1,nb_steps+1),'k-')
    end
    set(gca, 'YTick', unique([1, get(gca, 'YTick')]));
    set(gca, 'YTick', unique([2, get(gca, 'YTick')]));
    %legend('Radius','Level 1', 'Level 2')
    xlabel('Time','Fontweight','bold')
    ylabel('Radius of smallest enclosing circle','Fontweight','bold')
    title('Enclosing circle radius')
end

if memorise_pos_display && memorise_r_display
    saveas(gcf,['res/evolution/discrete_n_',num2str(n),'_square_size_',num2str(square_size),'_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'_seed_',num2str(seed)],'fig')
    saveas(gcf,['res/evolution/discrete_n_',num2str(n),'_square_size_',num2str(square_size),'_steps_',num2str(nb_steps),'_dt_',num2str(dt),'_v_',num2str(v),'_seed_',num2str(seed)],'epsc')
end



if save_conv_time_to_file
    for i=1:nb_trials
        fprintf(fileID, '%d\n', t_conv(i));
    end
end




if save_conv_time_to_file
   fclose(fileID);
end
