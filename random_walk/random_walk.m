

%%Setting up variables:
%n_particles = input("Number of Particles (insert number) = "); %number of particles
%timesteps = input("Timesteps (insert number) = "); %how many timesteps will we run for
%delta_r = input("Step size (insert number) = "); %set the distance it can move in a timestep
%position = zeros(2, n_particles, timesteps); %set up storage for position

n_particles = 10; %number of particles
timesteps = 200; %how many timesteps will we run for
delta_r = 1; %set the distance it can move in a timestep
position = zeros(2, n_particles, timesteps); %set up storage for position

%%Setting some general plot parameters:
sz = 55; %set markersize
%%Setting custom colors, my own preference. am into orange/blue combo lately
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";

%Saving plots as videos:
writerObj = VideoWriter(strcat("particles_randomwalk.avi")); %set the name of your video
writerObj.FrameRate = 25; %set the frame rate
open(writerObj); %open video writer

%Random walk; movement not constrained to a grid
figure
for i = 1:n_particles %loop over number of particles
     for j = 2:timesteps %loop over timesteps
            theta = 360 * rand(1, 1); %random angle between 0 and 2pi
            position(1, i, j) = position(1, i, j-1) + (delta_r) * cosd(theta); %x coord +x distance
            position(2, i, j) = position(2, i, j-1) + (delta_r) * sind(theta); %y coord +y distance
    end  
end 

%Plot 
for i = 1:timesteps %loop over timesteps
    %create a scatter plot for all particles over all timesteps
    sc = scatter(position(1, :, i), position(2, :, i), sz, "filled");
    %g = linspace(1,10,length(position(1, :, i)));
    %scatter(position(1, :, i), position(2, :, i), [], g, "filled");
    
    sc.MarkerEdgeColor = c_orange;
    sc.MarkerFaceColor = c_orange;
    sc.LineWidth = 0.5;
    
    %Custom legend:
    %params to display in the legend
	%legend_input = {strcat("{\it N} particles = " + ...
                        %num2str(n_particles) + "\newlineTimesteps = " + ...
                        %num2str(timesteps) + " au \newlineStep size = " + ...
                        %num2str(delta_r)) + " au"};
	%l = legend(legend_input); %show the legend
	%l.FontSize = 13 %legend font size
	%l.FontName = "Times" %legend fontface
	%legend('boxoff') %legend box on or off
   
	%Custom plot title:
    t = title({strcat("$N$ particles performing a random walk") ...
                      "(movement not constrained to a grid)"}); %title
	
	t.FontSize = 17; %plot title font size
    %t.FontColor = c_blue2;
    t.Interpreter = "latex";
    
    st = subtitle({strcat("{\it N} particles = " + ...
                          num2str(n_particles) + ", Timesteps = " + ...
                          num2str(timesteps) + ", Step size = " + ...
                          num2str(delta_r))});
                      
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
   
	% Other plot settings:
	xlim([-15 15]);
    ylim([-15 15]);
    
	box on %plot box on or off
	grid on %plot grid on or off
	axis square %axis aspect ratio (square or equal)
	pause(0.01) %increase or decrease the speed
   
	%Video settings:
	frame = getframe(gcf); %capture axes in a video  
	writeVideo(writerObj, frame); %write video
    %hold on %use hold on if you want to capture the particle trajectories
    %saveas(sc, sprintf('Q1_niko1_%d.svg', timesteps));
    
end 

 %Close video writer:
 close(writerObj);
 
