% This program is used to simulate the transmission of a fictional virus
% amongst N number of carriers in a confined space with set limits.
% WRITTEN BY: Imsara Samarasinghe
% DATE: 10/03/2020

% Housekeeping
clear
clc


%The simulation model has some fixed parameters chosen before hand
%---------------------MODEL PARAMETERS-------------------
Width = 1000;     % The width of the confined space along the x direction (m)
Height = 1000;    % The height of the confined space along the y direction (m)
Vmin = 0.1;       % The lower bound on the velocity of an individual in the simulation (m s^-1)
Vmax = 0.2;       % The upper bound on the velocity of an indiviual in the simulation (m s^-1)
T = 864000;       % The total time for which the simulation will run (s)
DeltaT = 10;      % The time step for the simulation (s)


%--------------------MODEL INPUTS------------------------
N = input('Please enter the number of individuals in the simulation:  ');


%--------------INITIAL INDIVIDUAL PROPERTIES-------------
position = rand(N,2).*(Width);                                                    % Initial x and y coordinates of the subjects    
VelocityMagnitudes = 0.1 + (0.2-0.1) .* rand(N,1);                                % Magnitude of velocity of each subject
theta = rand(N,1).*2.*pi;                                                         % Initial heading of each subject
VelocityComp = [VelocityMagnitudes.*cos(theta) VelocityMagnitudes.*sin(theta)];   % Components of velocity along the x and y directions 
   
% Finding the patient zero and assigning to an array
infected = zeros(N,1);     % Initialising the array of infected individuals
unlucky = randi([0 N],1);  % Finding the unlucky patient zero that is initially infected by the virus
infected(unlucky) = 1;     % Assigning that in the array infected

% Finding the positons of the healthy and unhealthy subjects
infected_Pos = zeros(1,2);
healthy_Pos = zeros(N,2);
for i = 1 : N
    if infected(i) == 1
        infected_Pos(:,1:2) = position(i,1:2);
    else
        healthy_Pos(i,1:2) = position(i,1:2);
    end
end
healthy_Pos( all(~healthy_Pos,2), : ) = []; % To get rid of zero row created during the loop.


% Plotting the initial locations of all the subjects
figure(1);
lineWidth = 2;
markerSize = 65;
scatter(infected_Pos(:,1),infected_Pos(:,2),markerSize,[1 0 0],'filled','MarkerEdgeColor','k'); hold on; 
scatter(healthy_Pos(:,1),healthy_Pos(:,2),markerSize,[0.1 0.6 0.1],'filled','MarkerEdgeColor','k');hold off; 
legend('Infected','Healthy','location','best');
legend('location','best');
xlim([0,Width]);
ylim([0,Width]);
axis square;
grid on;
box on;
hold off;
exportgraphics(figure(1),'InitialLocations.jpg'); % exporting the plot of the initial locations

%Setting up all arrays used within the loops
t=(linspace(0,T,T/DeltaT))./3600;    % Time vector
healthy = ~infected;         % Array of the healthy people
recovered = zeros(N,1);      % Array of the recovered people
incubation = 172800;         % Incubation period (s)
t_recover = 432000;          % Timetaken to recover (s)
t_infected = zeros(N,1);     % Array to keep track of time since being infected
sick = infected;             % Array of sick people who havent recovered yet    
asymptomatic = zeros(N,1);   % Array of infected but asymptomatic people
newPosition = zeros(N,2);    % New position array
asySum = zeros(N,1);         % sum of asymptomatic carriers
sickSum=zeros(N,1);          % sum of sick carriers
healSum=zeros(N,1);          % sum of healthy carriers
recSum=zeros(N,1);           % sum of recovered carriers

% For bar graph in the loop
X=categorical({'healthy','asymptomatic','sick','recovered'});
X=reordercats(X,{'healthy','asymptomatic','sick','recovered'});



for timeStep = 1 : (T/DeltaT) % The number of iterations the loop shall run for i.e. the main iteration
    newPosition = position+VelocityComp.*DeltaT; % new positions of the carriers
    for carrier1 = 1 : N
        if infected(carrier1)==1                                         %Only choosing the infected carriers
            if asymptomatic(carrier1)==1 && t_infected(carrier1)<incubation
                t_infected(carrier1)=t_infected(carrier1)+DeltaT;               %Time since being infected
            elseif asymptomatic(carrier1)==1 && t_infected(carrier1)==incubation
                asymptomatic(carrier1)=0;
                sick(carrier1)=1;
                t_infected(carrier1)=t_infected(carrier1)+DeltaT;
            elseif sick(carrier1)==1 && t_infected(carrier1)<t_recover
                t_infected(carrier1)=t_infected(carrier1)+DeltaT;
            elseif sick(carrier1)==1 && t_infected(carrier1)==t_recover
                infected(carrier1)=0; % Fully recovered so no longer infected
                sick(carrier1)=0;
                recovered(carrier1)=1;
            end
        end
        for carrier2 = 1:1:N
            if carrier2~=carrier1
                chance=rand; % Chance of infection for the carriers
                
                % Passing on the virus if carriers get within 2 meters of
                % healthy person
                if norm(newPosition(carrier1,:)-newPosition(carrier2,:))<=2 
                    if recovered(carrier2)||recovered(carrier1)         % if either recovered cannot pass on the virus
                        continue
                    elseif infected(carrier2)&&infected(carrier1)       % if both infected cannot pass on the virus
                        continue
                    elseif infected(carrier2)&& ~(infected(carrier1))
                        if sick(carrier2)
                            if chance<=0.5
                                infected(carrier1)=1;
                                asymptomatic(carrier1)=1;
                                healthy(carrier1)=0;
                            end
                        elseif asymptomatic(carrier2)
                            if chance<=0.3
                                infected(carrier1)=1;
                                asymptomatic(carrier1)=1;
                                healthy(carrier1)=0;
                            end
                        end
                    elseif infected(carrier1)&&~(infected(carrier2))
                        if sick(carrier1)
                            if chance<=0.5
                                infected(carrier2)=1;
                                asymptomatic(carrier2)=1;
                                healthy(carrier2)=0;
                            end
                        elseif asymptomatic(carrier1)
                            if chance<=0.3
                                infected(carrier2)=1;
                                asymptomatic(carrier2)=1;
                                healthy(carrier2)=0;
                            end
                        end
                    end
                end
            end
        end
        % Changing direction once at the walls using a user made function
        [newHeading,newComponents] = wallCollision(newPosition(carrier1,:),VelocityMagnitudes(carrier1,:),theta(carrier1));
        theta(carrier1)=newHeading;
        VelocityComp(carrier1,:)=newComponents;
        
    end
    % sum of all categories
    asySum(timeStep)=sum(asymptomatic);
    sickSum(timeStep)=sum(sick);
    recSum(timeStep)=sum(recovered);
    healSum(timeStep)=sum(healthy);

    currentTime=timeStep*DeltaT;   % Current time in seconds
    [Time]=dateTime(currentTime);  % Time vector output using a user made function
    
    % Structure
    sumCategory.healthy(timeStep)=sum(healthy);
    sumCategory.asymptomatic(timeStep)=sum(asymptomatic);
    sumCategory.sick(timeStep)=sum(sick);
    sumCategory.recovered(timeStep)=sum(recovered);
    sumCategory.time(timeStep,1:4)=Time(1,1:4);
    
       
    asyPos=zeros(N,2);   % Asymptomatic carriers locations
    sickPos=zeros(N,2);  % Sick carriers locations
    recPos=zeros(N,2);   % recovered carriers loactions
    healPos=zeros(N,2);  % healthy peoples locations
    for carrier1 = 1 : N
        % Splitting the positions according to the infection status
        if infected(carrier1)
            if asymptomatic(carrier1)
                asyPos(carrier1,1:2)=newPosition(carrier1,1:2);
            elseif sick(carrier1)
                sickPos(carrier1,1:2)=newPosition(carrier1,1:2);
            end
        elseif recovered(carrier1)
            recPos(carrier1,1:2)=newPosition(carrier1,1:2);
        elseif healthy(carrier1)
            healPos(carrier1,1:2)=newPosition(carrier1,1:2);
        end
    end
    % Removing all zero row entries created during above loop
    healPos( all(~healPos,2), : ) = [];
    recPos( all(~recPos,2), : ) = [];
    sickPos( all(~sickPos,2), : ) = [];
    asyPos( all(~asyPos,2), : ) = [];
    
    days=Time(1,1);
    hrs=Time(1,2);
    min=Time(1,3);
    sec=Time(1,4);
    b = mod(min,10);
    
    if b == 0
        figure(2)
        subplot(2,2,1)% Location
        markerSize=50;
        % asymptomatic locations
        scatter(asyPos(:,1),asyPos(:,2),markerSize,[1 0.5 0],'filled','MarkerEdgeColor','k');hold on;
        % sick positions
        scatter(sickPos(:,1),sickPos(:,2),markerSize,[1 0 0],'filled','MarkerEdgeColor','k');hold on;
        % recovered positions
        scatter(recPos(:,1),recPos(:,2),markerSize,[0 0 1],'filled','MarkerEdgeColor','k');hold on;
        % healthy Positions
        scatter(healPos(:,1),healPos(:,2),markerSize,[0.1 0.6 0.1],'filled','MarkerEdgeColor','k');hold off;  
        xlim([0 Width]);
        ylim([0 Height]);
        axis square
        grid on
        box on
        titleString=strcat('Days= ',num2str(Time(1)),' Hours= ',num2str(Time(2)));
        title(titleString);
        
        subplot(2,2,2) % Bar plot
        Y=[sum(healthy) sum(asymptomatic) sum(sick) sum(recovered)];
        b=bar(X,Y);
        b.FaceColor='flat';
        b.CData(1,:)=[0.1 0.6 0.1];b.CData(2,:)=[1 0.5 0];b.CData(3,:)=[1 0 0];b.CData(4,:)=[0 0 1];
        title(titleString);
        
        if hrs==0 && min==0 && sec==0 % history plot updates only after a day to keep the speeds up
            subplot(2,2,[3,4])% History plot
            plot(t(1:timeStep),healSum(1:timeStep),'Color',[0.1 0.6 0.1],'LineWidth',2);hold on;
            plot(t(1:timeStep),asySum(1:timeStep),'Color',[1 0.5 0],'LineWidth',2);hold on;
            plot(t(1:timeStep),sickSum(1:timeStep),'r','LineWidth',2);hold on;
            plot(t(1:timeStep),recSum(1:timeStep),'b','LineWidth',2);hold on;
            legend('healthy','asymptomatic','sick','recovered')
            xlabel('hours')
            ylabel('total')
            xlim([0 240])
            ylim([0 N])
        end
        
        
        % exporting plots at 2,4 and 6 days
        if days==2 && hrs==0 && min==0 && sec==0
            exportgraphics(figure(2),'2_Days.jpg')
        elseif days==4 && hrs==0 && min==0 && sec==0
            exportgraphics(figure(2),'4_Days.jpg')
        elseif days==6 && hrs==0 && min==0 && sec==0
            exportgraphics(figure(2),'6_Days.jpg')
        end

    end
    position=newPosition;      
end
% Daily summary
for j=1:86400
    if sumCategory.time(j,2)==0 && sumCategory.time(j,3)==0 && sumCategory.time(j,4)==0
        diary on
        diary dailySummary.txt
        disp(strcat('Day Number ',num2str(sumCategory.time(j,1))));
        disp(strcat('The number of healthy people= ',num2str(sumCategory.healthy(j))));
        disp(strcat('The number of sick people= ',num2str(sumCategory.sick(j))));
        disp(strcat('The number of asymptomatic people= ',num2str(sumCategory.asymptomatic(j))));
        disp(strcat('The number of recovered people= ',num2str(sumCategory.recovered(j))));
        diary off
    end      
end