clear; clc;
format long g;
char brd_fname;
% The following script reads ephemeris data and calculates the required kepler elements
% clock corrections etc. in order to supply the user with satellite navigation data
% such as x, y, z coordinates for all the available satellites given a starting time(epoch)
% up until the end of the RINEX file.
% Written by locateTempUserBash
% visit the user page @ github for further information
% or email using locateTempUserBash@yahoo.com

%------ test 08.06.2017 ------
    

brd_orb = Coord_Eph();
brd_fname='auto0010.05n';
date_in=[2005,1,1,0,0,0];
epoch=date_time(date_in);
time_interval = 15*60*60; % time interval between each epoch
brd_orb = broadcast_position('auto0010.05n',time_interval,brd_orb);
