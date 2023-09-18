%% PHOTOELASTICITY EXPERIMENT: isoclinic simulation

% This MATLAB file "isoclinic_simulation_code.m" contains code for simulating 
% isoclinic patterns and animating isoclinic changes on a photoelastic disc, 
% according to two theoretical models of physical stress maps - the Ovals model 
% and the Wulff net model, as explained in the "photoelasticity_report.pdf" file 
% (sections 1,2,4) and its appendix presentation file "photoelasticity report - 
% video appendix.ppsx".
% The "Photoelasticity presentation.ppsx" contains further demonstrations and 
% animations of the isoclinic simlulation process.
%%

% Code, paper and presentations were created as part of a photoelsticity experiment, 
% carried out in LAB B course of our undergarduate physics studies, at the Hebrew 
% University of Jerusalem.
%
% Created by Shlomo Danziger and Guy Elfersy (c) 
% 
% Comments, questions, etc. are welcome:
% ShlomoMD@gmail.com
%%

%% ADJUSTABLE GENERAL SETTINGS

clear; % clear workspace upon running (recommended). may be disabled

task = 1; 
    % choose TASK to be carried out:
    % (0): create isoclinic map, (1): display isoclinic change animation
model = 1; 
    % choose MODEL of disc stress map: (0) Ovals model, (1) Wulff net model
lines = 1; 
    % set to 1 to display isoclinic LINES, 0 otherwise
arrows = 1; 
    % set to 1 to display isoclinic ARROWS, 0 otherwise
polarizer = 0; 
    % set to 1 to display turing POLARIZER in background (recommended for 
    % "display isoclinic change animation" task), 0 otherwise
continuous = 1; 
    % set to 1 to loop animations CONTINUOUSLY (recommended for "display 
    % isoclinic change animation" task), 0 otherwise
save_pics = 0;
    % Set to 1 to SAVE PIC of plot in each iteration as a separate file in 
    % local folder (to create a GIF file independently from all pics later)
create_gif = 0;
    % set to 1 to automatically create a GIF file from all plots (the "gif" 
    % function must first be downloaded from File Exchange, for this option)
display_axes = 0;
    % set to 0 to get rid of axes (there's not much need for them in this 
    % simulation), 1 to display axes
    
%% MORE ADJUSABLE (AND NON-ADJUSATABLE) PARAMETERS

% ADJUSTABLE parameters for ARROWS:
max_head_size = 0.05;
auto_scale_factor = 1; % val recommended for isoclinic change animation: 1
length_factor = 0.1; % val recommended for isoclinic change animation: 0.1
line_width = 0.5;

% ADJUSTABLE parameters for ANGLES (corresponding to isoclinics):
number_of_angles = 18; 
    % number of isoclinics to appear in each quadrant
start_angle = 3; 
    % in degrees. may be changed according to angle of first sample taken 
    % at lab (note that angles 0,90 will cause division by 0)

% NON-ADJUSTABLE parameters for angles (corresponding to isoclinics):
angle_jump = 90/number_of_angles;
end_angle = start_angle + 90 - angle_jump;
angle = horzcat(start_angle:angle_jump:end_angle);
t_q1 = tan(angle*pi/180);
t_q2 = tan((angle+90)*pi/180);

% ADJUSTABLE PAUSE interval before appearance of next isoclinic:
pause_interval = 0.05;
    % make this small if number of angles is large

% ADJUSTABLE parameters for WULFF NET model:
n = 10;
    % higher number creates more (relatively) wide circles => shorter arrows
    
% ADJUSTABLE parameters for OVALS model:
num_of_ovals_determiner = 30; 
    % higher number creates more narrow ovals
a_jump = 2; 
    % lower number creates more (relatively) wide ovals => shorter arrows

% NON-ADJUSTABLE parameters for Ovals model:
a=0:a_jump/10:0.5*(num_of_ovals_determiner-1);
m=2.^a;

% ADJUSTABLE parameters for DISC
r = 10; %radius  (if axes are turned off then r's value doesn't matter)

% NON-ADJUSTABLE parameters for POLARIZER symbol
R_pol_in=r*1.2; % internal radius
R_pol_out=r*1.4;  % external radius

% COLORS
c_animation = 'r'; % color of isoclinics in isoclinic change animation
c_all = hsv(number_of_angles); % create color pallete for isoclinic map

%% MAIN

while true % if continuous mode is enabled, animation will run until CTRL+C is pressed
    
    hold off
    draw_disc(r);
    
    for i = 1:1:number_of_angles

        % draw turning polarizer in background (if enabled)
        if polarizer
            if i == 1 % initiate parameters for draw_polarizer(...) ahead
                l1=0; l2=0; l3=0; l4=0;
            end 
            [l1,l2,l3,l4] = draw_polarizer(angle,i,R_pol_in,R_pol_out,l1,l2,l3,l4);        
        end
        
        % set color(s) and clear screen if necessary, according to task
        if task == 0   % "create isoclinic map" task
            c = c_all(i,:); % color spectrum
        else           % "display isoclinic change animation" task
            hold off;
            draw_disc(r);
            c = c_animation;
        end 
        
        if display_axes == 0  % do not display axes unless requested
            axis off
        end
        
        % get coordinates of isoclinics, depending on model
        if model    % Wulff net model chosen
            [x_neg,y_of_x_neg,x_pos,y_of_x_pos] = get_coor_of_ith_isoclinic_Wulff_net(t_q1,t_q2,i,r,n);
        else        % Ovals model chosen
            [x_neg,y_of_x_neg,x_pos,y_of_x_pos] = get_coor_of_ith_isoclinic_ovals(t_q1,t_q2,i,r,m);
        end
    
        % plot isoclinic line (if enabled)
        if lines
            plot_isoclinic_line(x_neg,y_of_x_neg,x_pos,y_of_x_pos,c);
        end
    
        % plot isoclinic arrows (if enabled)
        if arrows
            plot_isoclinic_arrows(angle,i,x_neg,y_of_x_neg,x_pos,y_of_x_pos,c,r,polarizer,max_head_size,auto_scale_factor,length_factor,line_width);
        end

        % create graph, according to task
        if task == 0    % "create isoclinic map" task
            colormap(c_all);
            ylabel(colorbar, 'Angle (degrees) of disc relative to polarizers');
            caxis([0 90]);
        else            % "display isoclinic change animation" task
            title(['Angle (degrees) of disc relative to polarizers: ' num2str(angle(i))])
        end 

        if save_pics 
            save_pics_in_separate_files(angle,i); 
        end
        
        if create_gif 
            create_gif_automatically(i);
        end

        pause(pause_interval); % pause before next isoclinic

    end  % end of i-th iteration of "FOR" loop

    if continuous == 0 
        break
    end

end


%% FUNCTIONS

function draw_disc(r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th);
    yunit = r * sin(th);
    plot(xunit, yunit,'b','Linewidth',1);  % color & thickness of circumference
    axis square
    hold on
end

function [new_l1,new_l2,new_l3,new_l4] = draw_polarizer(angle,i,R_pol_in,R_pol_out,l1,l2,l3,l4)
            % delete old polarizer lines
            if 1 < i
                delete(l1);
                delete(l2);
                delete(l3);
                delete(l4);
            end
            % draw new polarizer lines + re-draw polarizer
            x_in=R_pol_in*cosd(angle(i));
            x_out=R_pol_out*cosd(angle(i));
            y_in=R_pol_in*sind(angle(i));
            y_out=R_pol_out*sind(angle(i));
            new_l1 = plot([x_out x_in],[y_out y_in],'black','Linewidth',0.5);
            new_l2 = plot([-x_out -x_in],[-y_out -y_in],'black','Linewidth',0.5);
            new_l3 = plot([y_out y_in],[-x_out -x_in],'black','Linewidth',0.5);
            new_l4 = plot([-y_out -y_in],[x_out x_in],'black','Linewidth',0.5);
            plot(R_pol_out*cos(0:pi/50:2*pi), R_pol_out*sin(0:pi/50:2*pi), 'black','Linewidth',0.5);
        end

function [x_neg,y_of_x_neg,x_pos,y_of_x_pos] = get_coor_of_ith_isoclinic_ovals(t_q1,t_q2,i,r,m)
        x_neg = -(t_q1(i))*r./sqrt(m.^2+m*(t_q1(i))^2); 
        y_of_x_neg = sqrt(r.^2-m.*x_neg.^2);
            % coordinates of points on ovals in 2nd quad where slopes of 
            % ovals' tangents equal tan(angle(i)*pi/180). See explicit 
            % formula in separate pic file.
        x_pos = -(t_q2(i))*r./sqrt(m.^2+m*(t_q2(i))^2);
        y_of_x_pos = sqrt(r.^2-m.*x_pos.^2);
            % coordinates of points on ovals in 1st quad where slopes of 
            % ovals' tangents equal tan((angle(i)+90)*pi/180)
end

function [x_neg,y_of_x_neg,x_pos,y_of_x_pos] = get_coor_of_ith_isoclinic_Wulff_net(t_q1,t_q2,i,r,n)
        abs_m_neg_max = abs(r*t_q1(i));
        m_neg = 0:r/n:abs_m_neg_max;
        x_neg = m_neg-sqrt((r^2+m_neg.^2)*((t_q1(i))^2)/((t_q1(i))^2+1));
        y_of_x_neg = sqrt(r^2+2*m_neg.*x_neg-x_neg.^2);
             % coordinates of points on circle-arches in 2nd quad where 
             % slopes of circle' tangents equal tan(angle(i)*pi/180). See 
             % explicit formula in separate pic file.
        abs_m_pos_max = abs(r*t_q2(i));
        m_pos = 0:r/n:abs_m_pos_max;
        x_pos = -m_pos+sqrt((r^2+m_pos.^2)*((t_q2(i))^2)/((t_q2(i))^2+1));
        y_of_x_pos = sqrt(r^2-2*m_pos.*x_pos-x_pos.^2);
             % coordinates of points on circle-arches in 1st quad where 
             % slopes of circle' tangents equal tan((angle(i)+90)*pi/180)
end

function plot_isoclinic_line(x_neg,y_of_x_neg,x_pos,y_of_x_pos,color1)
            plot(x_neg,y_of_x_neg,'Linewidth',2,color=color1);
            plot(-x_neg,-y_of_x_neg,'Linewidth',2,color=color1);
            plot(x_pos,y_of_x_pos,'Linewidth',2,color=color1);
            plot(-x_pos,-y_of_x_pos,'Linewidth',2,color=color1);
end

function plot_isoclinic_arrows(angle,i,x_neg,y_of_x_neg,x_pos,y_of_x_pos,c,r,polarizer,max_head_size,auto_scale_factor,length_factor,line_width)
        if polarizer==0
            set(gca, 'XLim', [-1.05*r 1.05*r], 'YLim', [-r 1.1*r]); 
            % prevent graph from getting warped, if polarizer symbol is disabled
        end        
        u_neg = repelem(cosd(angle(i))*length_factor,length(x_neg));
        v_neg = repelem(sind(angle(i))*length_factor,length(y_of_x_neg));
        u_pos = repelem(cosd(angle(i)+90)*length_factor,length(x_pos));
        v_pos = repelem(sind(angle(i)+90)*length_factor,length(y_of_x_pos));
        h1 = quiver(x_neg-u_neg/2,y_of_x_neg-v_neg/2,u_neg,v_neg,auto_scale_factor, 'Linewidth',line_width,'MaxHeadSize',max_head_size,color=c);
        set(h1,'MaxHeadSize',max_head_size,'AutoScaleFactor',auto_scale_factor);
        h2 = quiver(-x_neg-u_neg/2,-y_of_x_neg-v_neg/2,u_neg,v_neg,auto_scale_factor, 'Linewidth',line_width,'MaxHeadSize',max_head_size,color=c);
        set(h2,'MaxHeadSize',max_head_size,'AutoScaleFactor',auto_scale_factor);
        h3 = quiver(x_pos-u_pos/2,y_of_x_pos-v_pos/2,u_pos,v_pos,auto_scale_factor, 'Linewidth',line_width,'MaxHeadSize',max_head_size,color=c);
        set(h3,'MaxHeadSize',max_head_size,'AutoScaleFactor',auto_scale_factor);
        h4 = quiver(-x_pos-u_pos/2,-y_of_x_pos-v_pos/2,u_pos,v_pos,auto_scale_factor, 'Linewidth',line_width,'MaxHeadSize',max_head_size,color=c);
        set(h4,'MaxHeadSize',max_head_size,'AutoScaleFactor',auto_scale_factor); 
end

function save_pics_in_separate_files(angle,i)
        saveas(gcf,['Isoclinic_map_after_angle(deg)_' num2str(angle(i)) '.jpg']);
end

function create_gif_automatically(i)
        if i == 1
            gif('Isoclinic_map_creation.gif','DelayTime',0.001)
        else
            gif
        end
end
