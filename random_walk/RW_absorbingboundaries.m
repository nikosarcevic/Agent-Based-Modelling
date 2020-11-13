%%PHY8005 | Module 1 | Agent Based Modelling 
%%Workshop 1 | Question 1-2
%%Nikolina Sarcevic
%%n.sarcevic2@newcastle.ac.uk
%%student number 200611321
%%-------------------------------------------------------------------------

%Setting up variables:
%n_particles = input("Number of Particles (insert number) = "); %number of particles
%timesteps = input("Timesteps (insert number) = "); %how many timesteps will we run for
%delta_r = input("Step size (insert number) = "); %set the timestep displacement
%abs_bound = input("Absorbing boundary value (insert number) = "); %set absorbing boundary
%position = zeros(2, n_particles, timesteps); %set up storage for position

n_particles = 10; %number of particles
timesteps = 100; %how many timesteps will we run for
delta_r = 1; %set the distance it can move in a timestep
position = zeros(2, n_particles, timesteps); %set up storage for position
abs_bound = 5; %set absorbing boundary

%Setting some general plot parameters:
sz = 55; %set markersize
%setting custom colors, my own preference. am into orange/blue combo lately
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";
%position of the box with appropriate absorbing boundaries
rectange_pos = [-abs_bound -abs_bound abs_bound*2 abs_bound*2];

%Saving plots as videos:
writerObj = VideoWriter(strcat("particles_abs_bounds.avi")); %set the name of your video
writerObj.FrameRate = 25; %set the frame rate
open(writerObj); %open video writer

%Random walk; movement not constrained to a grid; absorbing boundaries included
figure
for i=1:n_particles %loop over particles
    for j=2:timesteps %loop over timesteps
        %boundary condition for x and y coord greater than set boundary
        if abs(position(1, i, j-1)) > abs_bound || abs(position(2, i, j-1)) > abs_bound
            position(1, i, j) = position(1, i, j-1); %for x coord
            position(2, i, j) = position(2, i, j-1); %for y coord
        else
            theta=360*rand(1, 1); %random angle between 0 and 2pi
            position(1, i, j) = position(1, i, j-1) + (delta_r) * cosd(theta); %x coord +delta_x distance  
            position(2, i, j) = position(2, i, j-1) + (delta_r) * sind(theta); %y coord +delta_y distance
        end
    end  
end 

%Plot
for i=1:timesteps %loop over timesteps
    %create a scatter plot for all particles over all timesteps
	sc = scatter(position(1, :, i),position(2, :, i), sz);
    sc.MarkerEdgeColor = c_orange;
    sc.MarkerFaceColor = c_orange;
    sc.LineWidth = 0.5;
        
	%Custom legend:
    %params to display in the legend
	%legend_input = {strcat("{\it N} particles = " + ...
                    %num2str(n_particles) + "\newlineTimesteps = " + ...
                    %num2str(timesteps) + " au \newlineStep size = " + ...
                    %num2str(delta_r) + " au \newlineBoundary length = " + ...
                    %num2str(2*abs_bound) + " au")};
	%l = legend(legend_input); %show the legend
	%l.FontSize = 13 %legend font size
	%l.FontName = "Times" %legend fontface
	%legend('boxoff') %legend box on or off
   
	%Custom plot title:
	t = title({"$N$ particles performing a random walk",...
        "(with absorbing boundary condition)"}); %title
    
	t.FontSize = 17; %plot title font size
    t.FontName = "Times";
    %t.Color = c_blue2;
    t.Interpreter = "latex";
    
    st = subtitle({strcat("{\it N} particles = " + ...
                          num2str(n_particles) + ", Timesteps = " + ...
                          num2str(timesteps) + ", Step size = " + num2str(delta_r)) + "," ...
                          ("Absorbing boundary ($(x_1, x_2) (y_1, y_2)$) =" + ...
                          "((" + num2str(-abs_bound) + "," + num2str(abs_bound) + ")" + ...
                          "(" + num2str(-abs_bound) + "," + num2str(abs_bound) + "))")});
                      
    st.FontName = "Times";
    st.FontSize = 15;
    st.Interpreter = "latex";
    %st.FontColor = c_blue1;
        
	%x and y label settings:
	x = xlabel("\fontname{Times}{\it x} coordinate [au]"); %x label title
	y = ylabel("\fontname{Times}{\it y} coordinate [au]"); %y label title
	x.FontSize = 15; %x label font size
	y.FontSize = 15; %y label font size
    x.FontName = "Times" %x fontface
    y.FontName = "Times" %y fontface
    
    %Absorbing boundary indicator:
    %plot a rectangle that corresponds to our abosrbing boundaries
    rec = rectangle("Position", rectange_pos,...
              "EdgeColor", c_blue1,...
              "LineWidth",1.5,...
              "LineStyle", "-");
	
	%Other plot settings:
	xlim([-abs_bound-10 abs_bound+10]); ylim([-abs_bound-10 abs_bound+10]); %plot lims
	box on %plot box on or off
	grid on %plot grid on or off
	axis square %axis aspect ratio (square or equal)
	pause(0.01) %increase or decrease the speed
        
    %Video settings:
    frame = getframe(gcf); %capture axes in a video  
    writeVideo(writerObj, frame); %write video
    %hold on %use hold on if you want to capture the particle trajectories
    %saveas(sc, sprintf('Q2_niko1_%d.svg', timesteps));
end

 %Close video writer:
 close(writerObj);
    