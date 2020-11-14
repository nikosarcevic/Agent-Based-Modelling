

%Setting up variables:
%n_particles = input("Number of Particles (insert number) = "); %number of particles
%timesteps = input("Timesteps (insert number) = "); %how many timesteps will we run for
%delta_r = input("Step size (insert number) = "); %set the timestep displacement
%per_bound = input("Periodic boundary value (insert number) = "); %set periodic boundary
%position = zeros(2, n_particles, timesteps); %set storage for position

n_particles = 10; %number of particles
timesteps = 100; %how many timesteps will we run for
delta_r = 1; %set the timestep displacement
per_bound = 5; %set periodic boundary
position = zeros(2, n_particles, timesteps); %set storage for position

%Setting some general plot parameters:
sz = 55; %set markersize
%setting custom colors, my own preference. am into orange/blue combo lately
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";
%position of the box with appropriate absorbing boundaries
%rectangle_pos = [-per_bound -per_bound per_bound*2 per_bound*2];

%Saving plots as videos:
writerObj = VideoWriter(strcat("particles_per_bounds.avi")); %set the name of your video
writerObj.FrameRate = 25; %set the frame rate
open(writerObj); %open video writer
rectangle_pos = [-per_bound -per_bound per_bound*2 per_bound*2]; %position of the box with appropriate periodic boundaries

%Random walk; movement not constrained to a grid; periodic boundaries included
figure
for i = 1:n_particles
    for j = 2:timesteps
        theta = 360*rand(1, 1); %random angle between 0 and 2pi
        position(1, i, j) = position(1, i, j-1) + (delta_r) * cosd(theta); %x coord +x displacement  
        position(2, i, j) = position(2, i, j-1) + (delta_r) * sind(theta); %y coord +y displacement
            
        %conditions when particle crosses the periodic boundary
        if abs(position(1, i, j)) > per_bound %x coord when x>boundary
            if position(1, i, j) < 0 %less than zero case
               position(1, i, j) = position(1, i, j) + 2 * per_bound;
             elseif position(1, i, j) > 0 %greater than zero case
                    position(1, i, j) = position(1, i, j) - 2 * per_bound;
                end
             elseif abs(position(2, i,j )) > per_bound %y coord when y>boundary
                 if position(2, i, j) < 0 %less than zero case
                    position(2, i, j) = position(2, i, j) + 2 * per_bound;
                elseif position(2, i, j) > 0 %greater than zero case
                    position(2, i, j) = position(2, i, j) - 2 * per_bound;
                end
            else
            end
        end    
end

%Plot
for i=1:timesteps %loop over timesteps
    %create a scatter plot for all particles over all timesteps
	sc = scatter(position(1, :, i),position(2, :, i), sz,...
            "filled",'MarkerEdgeColor', c_orange,...
            "MarkerFaceColor", c_orange,...
            "LineWidth", 0.5)
        
	%Custom legend:
    %params to display in the legend
	%legend_input = {strcat("{\it N} particles = " + ...
                    %num2str(n_particles) + "\newlineTimesteps = " + ...
                    %num2str(timesteps) + " \newlineStep size = " + ...
                    %num2str(delta_r) + " \newlineBoundary length = " + ...
                    %num2str(2*per_bound))};
	%l = legend(legend_input); %show the legend
	%l.FontSize = 13 %legend font size
	%l.FontName = "Times" %legend fontface
	%legend('boxoff') %legend box on or off
    
    %Custom plot title:
	t = title({"$N$ particles performing a random walk",...
        "(with periodic boundary condition)"}); %title
    
	t.FontSize = 17; %plot title font size
    t.FontName = "Times";
    %t.Color = c_blue2;
    t.Interpreter = "latex";
    
    st = subtitle({strcat("{\it N} particles = " + ...
                          num2str(n_particles) + ", Timesteps = " + ...
                          num2str(timesteps) + ", Step size = " + num2str(delta_r)) + "," ...
                          ("Periodic boundary ($(x_1, x_2) (y_1, y_2)$) =" + ...
                          "((" + num2str(-per_bound) + "," + num2str(per_bound) + ")" + ...
                          "(" + num2str(-per_bound) + "," + num2str(per_bound) + "))")});
                      
    st.FontName = "Times";
    st.FontSize = 15;
    st.Interpreter = "latex";
    %st.FontColor = c_blue1;
        
	%x and y label settings:
	x = xlabel("\fontname{Times}{\it x} coordinate"); %x label title
	y = ylabel("\fontname{Times}{\it y} coordinate"); %y label title
	x.FontSize = 15; %x label font size
	y.FontSize = 15; %y label font size
    x.FontName = "Times" %x fontface
    y.FontName = "Times" %y fontface
    
    %Periodic boundary indicator:
    %plot a rectangle that corresponds to our periodic boundaries
    rectangle("Position", rectangle_pos,...
              "EdgeColor", c_blue1,...
              "LineWidth",1.5,...
              "LineStyle", "-");
	
	%Other plot settings:
	xlim([-per_bound-10 per_bound+10]); %x lims
    ylim([-per_bound-10 per_bound+10]); %y lims
	box on %plot box on or off
	grid on %plot grid on or off
	axis square %axis aspect ratio (square or equal)
	pause(0.01) %increase or decrease the speed
        
    %Video settings:
    frame = getframe(gcf); %capture axes in a video  
    writeVideo(writerObj, frame); %write video
    %hold on %use hold on if you want to capture the particle trajectories
   
end
    
 %Close video writer:
 close(writerObj);
    
