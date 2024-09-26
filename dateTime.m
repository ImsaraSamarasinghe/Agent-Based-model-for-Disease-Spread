function [Time]=dateTime(sec)
% DESCRIPTION
% Accepts the time in seconds and outputs in the days,hours,minutes,
% seconds format.
%
% INPUTS
% sec= The time in seconds
%
% OUTPUTS
% Time= Array in Days Hours Minutes Seconds format

day=fix(sec/(24*3600));

sec= mod(sec,(24*3600));
hour=fix(sec/3600);

sec=rem(sec,3600);
minutes=fix(sec/60);

sec=rem(sec,60);
seconds=sec;
Time=[day hour minutes seconds];


