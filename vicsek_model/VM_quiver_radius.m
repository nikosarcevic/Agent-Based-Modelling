%%PHY8005 | Module 1 | Agent Based Modelling 
%%Workshop 2 | Question 2-2
%%Nikolina Sarcevic
%%n.sarcevic2@newcastle.ac.uk
%%student number 200611321
%%-------------------------------------------------------------------------

% Setting up simulation parameters:
n_particles = 200; % number of particles
t_max = 10; % duration
dt = 1; % time step
v = 0.5; % absolute velocity
L = 10; % cell size
noise = 0.0; % noise initial value)
r = 0.0; % interacting radius (initial value)
per_bound = 1; % 1: periodic boundary condition, 0: unlimited
rho = n_particles/L^2 % density
n_sims = 10; %number of simulations

%n_sims = input('Enter number of simulations:');
v_a = zeros(2, n_sims); % an empty array for the absolute velocity v_a

% Setting custom colors
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";

% Saving plots as videos: 
%writerObj = VideoWriter(strcat("vicsek_quiver_noise_const.avi")); % set the name of your video
%writerObj.FrameRate = 25; % set the frame rate
%open(writerObj); % open video writer

for i = 1:n_sims % loop over the number of simulations
    
    noise = 1; % set interaction radius to a constant value
	v_a(1, i) = r; % get noise values
	v_a(2, i) = vicsek(n_particles, dt, t_max, v, per_bound, L, r, noise); % get alignment values
    r = r + 0.1; % increase the the interacting radius by .1 in each iteration
    
end

% Define a Vicsek model and return an absolute velocity v_a
function [v_a] = vicsek(n_particles, dt, t_max, v, per_bound, L, r, noise)

x = L * rand(1, n_particles); % random initial position of x coord
y = L * rand(1, n_particles); % random initial position of y coord
theta = 2 * pi * (rand(1, n_particles) - 0.5); % random initial direction angle
rho = n_particles/L^2 % density

% Custom colors
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";
c_gray1 = "#808080";
c_gray2 = "#DCDCDC";

    for time = 1:t_max % loop from t_initial to t_max
            D = pdist([x' y'], 'euclidean'); % return the Euclidean distance between pairs of observations

            if per_bound == 1 % periodic boundary conditions:
                % boundary conditions for x coord
                tmp_x(x < r) = L + x(x < r);                    
                tmp_x(x > L-r) = x(x > L-r) - L;                      
                tmp_x(r <= x & x <= L-r) = x(r <=x & x <= L-r);       
                
                % boundary conditions for y coord
                tmp_y(y < r) = L + y(y < r);                       
                tmp_y(y > L-r) = y(y > L-r)-L;                     
                tmp_y(r <=y & y <= L-r) = y(r <= y & y <= L-r);       

                tmp_D = pdist([tmp_x' tmp_y'], 'euclidean'); % pairwise distance between pairs of observations      
                D = min([D; tmp_D]); % get minimum elements of an array
            end

            M = squareform(D); % format distance matrix
            [l1, l2]=find(0 < M & M < r); % find indices and values of nonzero element

            % calculate a circular mean for theta
            for i = 1:n_particles
                list = l1(l2 == i);
                if ~isempty(list)
                    ave_theta(i) = atan2(mean(sin(theta(list))), mean(cos(theta(list))));
                else
                    ave_theta(i) = theta(i);
                end
            end
        
            % get  the four quiver plot components
            x = x + v * cos(theta) * dt;
            y = y + v * sin(theta) * dt;
            U = v * cos(theta) * dt;
            V = v * sin(theta) * dt; 

            % periodic boundary conditions for the box:
            if per_bound == 1
                x(x < 0) = L + x(x < 0);
                x(L < x) = x(L < x) - L;
                y(y < 0) = L + y(y < 0);
                y(L < y) = y(L < y) - L;
            end                              
            
            % get velocity values
            v_comp = zeros(n_particles, 2); % empty array for velocity components v_x and v_y for all particles
            v_comp(:, 1) = v * cos(theta); % store v_x component for all particles
            v_comp(:, 2) = v * sin(theta); % store v_y component for all particles
            v_xy = sum(v_comp); % sum of all velocity components
            v_mag = sqrt(v_xy(1)^2 + v_xy(2)^2); % velocity magnitude
            
            theta = ave_theta + noise*(rand(1, n_particles) - 0.5); % angle of movement with the noise term correction
            v_a = (v_mag/(n_particles * v)); % absolute velocity (avg of all velocity vectors of all agents)
        
            % Quiver plot:
            q = quiver(x, y, U, V); 
            q.Color = c_orange % arrow color
            
            xlim([0 L]); % plot limits (0 < x < L)
            ylim([0 L]); % plot limits (0 < y < L)

            % Custom plot title settings:
            ptitle = title({strcat("Viscek model for " + ...
                                   num2str(n_particles) + " particles") ...
                                   ("with an increasing interacting radius") ...
                                   "and a periodic boundary condition"}); % title with preferred information
                    
            ptitle.FontSize = 17; % plot title font size
            %ptitle.Color = c_blue1; % plot title font color
            ptitle.FontName = "Times"; % plot title fontface
            ptitle.Interpreter = "latex" % plot title interpreter

            % Custom plot subtitle settings:
            stitle = subtitle({strcat("$\rho$ = " + num2str(rho) + ...
                                      ", $r_i$ = " + num2str(r/t_max) + ...
                                      ", $\eta$ = " + num2str(noise) + ...
                                      ", $\mathrm{d}t$ = " + num2str(dt) + ...
                                      ", $t_\mathrm{max}$ = " + num2str(t_max))}); % subtitle with preferred information params

            stitle.FontSize = 15; % plot subtitle font size
            stitle.FontName = "Times"; % plot subtitle fontface
            %stitle.Color = c_blue1; % plot subtitle color
            stitle.Interpreter = "latex"; % interpreter

            % x and y label settings:
            x_ax = xlabel("\fontname{Times}{\it x} coordinate"); %x label title
            y_ax = ylabel("\fontname{Times}{\it y} coordinate"); %y label title
            x_ax.FontSize = 15; %x label font size
            y_ax.FontSize = 15; %y label font size
            x_ax.FontName = "Times" %x fontface
            y_ax.FontName = "Times" %y fontface
            
            axis square % set axis to square
            pause(0.01) % speed
                saveas(q, sprintf("VM_quiver_radius%d.svg", time)) % save quiver plot figs as svg for each time iteration
    end

end

