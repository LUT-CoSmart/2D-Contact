function [Fc,Kc,GapNab,Gap] = Contact(Body1,Body2,approach,ContactPointfunc,Gapfunc)

% Initialize the gap 
Gap = 0;

% Initialize the global contact forces
Fc = zeros(Body1.nx + Body2.nx,1);

% Initialize the stiffnesses
Kc = zeros(Body1.nx + Body2.nx, Body1.nx + Body2.nx);

% Initialize the gap variation
GapNab = zeros(Body1.nx + Body2.nx,1);

% Initialize the gap distribution 
% GapDOFs = zeros(Body1.nx + Body2.nx,1);

if approach.Type ~= "None" % we have contact algorithm    
    addpath("Contact\ProjectionFunctions")
    [Fc,Kc,GapNab,Gap] = ContactVariation(Body1,Body2,approach,ContactPointfunc,Gapfunc);    
    
end


