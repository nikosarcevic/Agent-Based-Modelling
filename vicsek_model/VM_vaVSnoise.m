
% Setting up simulation parameters:
n_particles = 200; % number of particles
t_max = 100; % duration
dt = 1; % time step
v = 0.5; % initial velocity
L = 10; % cell size
noise = 0.0; % noise initial value)
r = 0.0; % interacting radius (initial value)
per_bound = 1; % 1: periodic boundary condition, 0: unlimited
n_sims = 100; % number of simulations
rho = n_particles/L^2 % density

%n_sims = input('Enter number of simulations:');
v_a = zeros(2, n_sims); % an empty array for the alignment angle v_a

% Setting custom colors, my own preference. am into orange/blue combo lately
c_blue1 = "#03396c";
c_blue2 = "#005b96";
c_blue3 = "#6497b1";
c_orange = "#ff7400";
c_gray1 = "#808080";
c_gray2 = "#DCDCDC";

for i = 1:n_sims % loop over the number of simulations
    r = 1; % set interaction radius to a constant value
	v_a(1, i) = noise; % get noise values
	v_a(2, i) = vicsek(n_particles, dt, t_max, v, per_bound, L, r, noise); % get absolute velocity values (instanteneous order param)
    noise = noise + 0.1; % increase the noise by .1 in each iteration
    
end

% Plot the aligment angle v_a as a function of the increasing noise:
plt = plot(v_a(1, :),v_a(2, :)); %, "-gs"

plt.Color = c_orange; % line color
plt.LineWidth = 1.5; % line width
plt.MarkerEdgeColor = c_blue1; % marker color

% Custom plot title settings:
ptitle = title({strcat("Viscek model for " + ...
                    num2str(n_particles) + " particles") ...
                    ("with an increasing noise term") ...
                    "and a periodic boundary condition"}); % title with preferred information
                                
ptitle.FontSize = 17; % plot title font size
ptitle.FontName = "Times"; % plot title fontface
%ptitle.Color = c_gray1;% plot title font color
ptitle.Interpreter = "latex"; % interpreter

% Custom plot subtitle settings:
stitle = subtitle({strcat("$\rho$ = " + num2str(rho) + ...
                          ", $r$ = " + num2str(r) + ...
                          ", $\eta_i$ = " + num2str(noise/t_max) + ...
                          ", $\mathrm{d}t$ = " + num2str(dt) + ...
                          ", $t_\mathrm{max}$ = " + num2str(t_max))}); % subtitle with preferred information params
                      
stitle.FontSize = 15; % plot subtitle font size
stitle.FontName = "Times"; % plot subtitle fontface
%stitle.Color = c_gray1; % plot subtitle color
stitle.Interpreter = "latex"; % interpreter

% Axes label settings:
x_ax = xlabel("Noise $\eta$"); % x label
y_ax = ylabel("Absolute velocity $v_\mathrm{a}$"); % y label

x_ax.FontSize = 15; % x label font size
y_ax.FontSize = 15; % y label font size
x_ax.FontName = "Times"; % x fontface
y_ax.FontName = "Times"; % y fontface
x_ax.Interpreter = "latex";
y_ax.Interpreter = "latex";

ylim([-0.1 1.1]); % set the limit for y axis, allow for some extra space (0 < v_a < 1)

% Add two horizontal lines for v_a min and max
yline1 = yline(0, "-", "not aligned"); % v_a min value
yline2 = yline(1, "-", "completely aligned"); % v_a max value

% Horizontal line attributes:
yline1.Color = c_blue1; 
yline2.Color = c_blue1;

yline1.LineWidth = 2;
yline2.LineWidth = 2;

yline1.FontName = "Times";
yline2.FontName = "Times";

yline1.FontSize = 15;
yline2.FontSize = 15;

yline1.LabelHorizontalAlignment = "left";

box on % include the plot box
grid on % plot the grid
axis square % keep the axes square
saveas(gcf,'VM_vaVSnoise.svg')

% Define a Vicsek model and return absolute velocity (an order parameter)
function [v_a] = vicsek(n_particles, dt, t_max, v, per_bound, L, r, noise)
x = L * rand(1, n_particles); % random initial position of x coord
y = L * rand(1, n_particles); % random initial position of y coord
theta = 2 * pi * (rand(1, n_particles) - 0.5); % random initial direction angle

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
            v_comp = zeros(n_particles, 2); % empty array for velocity components (v_x and v_y) for all particles
            v_comp(:, 1) = v * cos(theta); % store v_x component for all particles
            v_comp(:, 2) = v * sin(theta); % store v_y component for all particles
            v_xy = sum(v_comp); % sum of all velocity components
            v_mag = sqrt(v_xy(1)^2 + v_xy(2)^2); % velocity magnitude
            
            theta = ave_theta + noise*(rand(1, n_particles) - 0.5); % angle of movement with the noise term correction
            v_a = (v_mag/(n_particles * v)); % absolute velocity (avg of all velocity vectors of all agents)
        
    end

end
