function [newHeading,newComponents] = wallCollision(pos_new,VelocityMagnitudes,heading)
% DESCRIPTION
% This function accepts the position of a  subject within a domain and
% redirects them away from the walls of the domain if they come too close
% to the walls.
%
%INPUTS
%       pos_new - The position of each individual (An array containing
%                  both x and y coordinates).
%       VelocityMagnitudes - The velocity magnitudes of each individual
%OUTPUTS
%       Newheading - the new direction of travel assigned to the subject
%       newComponents - the new components of velocity applied to the
%                       individuals (an array containing u and v
%                       components)
newHeading=0;
if pos_new(1,1)<=2
    if (pi/2<=heading) && (heading<=1.5*pi)
        if rand<=0.5
            newHeading =(1.5+(2-1.5)*rand)*pi;
        else
            newHeading=0.5*rand*pi;
        end
    end
elseif pos_new(1,1)>=998
    newHeading=(0.5+(1.5-0.5)*rand)*pi;
    
elseif pos_new(1,2)<=2
    if (pi<=heading)&&(heading<=2*pi)
        newHeading=rand*pi;
    end
elseif pos_new(1,2)>=998
    if ~((pi<=heading)&&(heading<=2*pi))
        newHeading=(1+(2-1)*rand)*pi;
    end
else
    newHeading=heading;
end

newComponents=[VelocityMagnitudes*cos(newHeading) VelocityMagnitudes*sin(newHeading)];

        
