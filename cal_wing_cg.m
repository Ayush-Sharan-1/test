function [ x,y ] = cal_wing_cg(Croot,Ctip,LEangle,span,DFC )
%Calculates CG of the Wing wrt the LEMAC in feet
%   Detailed explanation goes here
LEangle=degtorad(LEangle);
y=(0.35*span)/2;
c=(2*y*(Ctip-Croot))/span+Croot
Xfs=0.25*c; %front spar locations, VARIABLE
Xrs=0.6*c; % rear spar location, VARIABLE
x=0.7*(Xrs-Xfs)+y*tan(LEangle)-DFC*tan(LEangle);
%git test
%git test 2
%git test 3

end

