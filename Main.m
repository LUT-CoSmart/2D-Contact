clc, clear, close all
format long;
addpath("ElementFunctions");
addpath("Forces")
addpath("Meshing")
addpath("ProcessAnalysis")
addpath("Contact")
%########## Reads element's data ###############################
ElementData;   
%########## Reads problem's data ###############################
ProblemData;
%########## Element positioning (from (0.0) coord. ) ###########
Body1.shift.x = 0;
Body1.shift.y = 0;

Body2.shift.x = 0;
Body2.shift.y = -Body2.Ly;

%#################### Mesh #########################################
dx1 = 16;
dy1 = 2;

dx2 = 8;
dy2 = 1;
%##################### Contact ############################
approach = 1; % 0 - none; 
              % 1 - penalty 
              % 2 - Nitsche (linear of gap), 3- Nitsche (nonlinear of gap), 4 - all items    
              % 5 - Lagrange multiplier    
              % 6 - penalty (simplified ): it's very simplified, even without gap redistribution over nodes              
              % 7 - Augumented Lagrange multiplier
              % 8 - Lagrange multiplier (nonlinear constrain)
              % 9 - perturbed Lagrangian method 
              % 10 - perturbed Lagrangian method (nonlinear constrain) 
% Hyperparameters 
pn = 1e10;  % penalty
PointsofInterest = "LinSpace"; % options: "nodes", "Gauss", "LinSpace" 
% N.B.: "LinSpace" with n == 2 is equal to "nodes"; 
% Number of "LinSpace" + 1 = number of n in "Gauss" ('cause the first point of elements is omitted)
n = 3; % number of points per segment (Gauss & LinSpace points)

ContactPointfunc  = ContactPointSetting(PointsofInterest,n);
Gapfunc = GapCalculationSetting(PointsofInterest, n);

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(DofsAtNode,Body1);
Body2 = CreateFEMesh(DofsAtNode,Body2);
% % %#################### BC  ###########################################
Body1.loc.x = 0; 
Body1.loc.y = 'all';  % Number (Location of nodes along the axis) or 'all' can be an option

Body2.loc.x = 0; 
Body2.loc.y = 'all'; 

Body1 = CreateBC(Body1);
Body2 = CreateBC(Body2); 

%##################### Loadings ######################
% local positions (assuming all bodies in (0,0) )
Body1.Fext.x = 0; 
Body1.Fext.y =-62.5*10^7;


Body1.Fext.loc.x = Body1.Lx;
Body1.Fext.loc.y = 'all';

Body2.Fext.y = 0; 
Body2.Fext.x = 0;

Body2.Fext.loc.x = Body2.Lx;
Body2.Fext.loc.y = 'all';

%##################### Egde nodes #########################
Body1.edge1.loc.x = Body1.Lx;
Body1.edge1.loc.y = 0;

Body1.edge2.loc.x = Body1.Lx;
Body1.edge2.loc.y = Body1.Ly;

Body2.edge1.loc.x = Body2.Lx;
Body2.edge1.loc.y = 0;

Body2.edge2.loc.x = Body2.Lx;
Body2.edge2.loc.y = Body2.Ly;
 
% Identification of possble contact surfaces
% local positions (assuming all bodies in (0,0) )
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;   
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;  
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);% 

%##################### Newton iter. parameters ######################
imax=40;
tol=1e-4;   
type = "cubic"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
steps= 20;

% %#################### Processing ######################
total_steps = 0;
titertot=0;  
for ii = 1:steps
    
        lambda_converged = false;
        if (approach == 5) || (approach == 8) || (approach == 9) || (approach == 10)
            lambda = 0;
        else
            lambda = zeros(Body1.ndof + Body2.ndof,1); % Lagrange item initiation
        end
       
        Body1 = CreateFext(ii,steps,Body1,type);
        Body2 = CreateFext(ii,steps,Body2,type);
   
        while (~lambda_converged) % special case for Augumented Lagrange            
            for jj = 1:imax % contact convergence
                tic;
                
                total_steps = total_steps + 1;
                
                % interacation of two bodies
                [Fc,Kc,GapNab,Gap] = Contact(Body1,Body2,pn,approach,ContactPointfunc,Gapfunc);
    
                % inner forces of the each body
                Body1 = Elastic(Body1);
                Body2 = Elastic(Body2);
                                       
                [uu_bc, deltaf, lambda] = Assemblance(Body1, Body2, Fc,Kc,GapNab,approach,pn,lambda);
                        
                % Displacement separation
                Body1.u(Body1.bc) = Body1.u(Body1.bc) + uu_bc(1:Body1.ndof);
                Body2.u(Body2.bc) = Body2.u(Body2.bc) + uu_bc(Body1.ndof + 1:Body1.ndof + Body2.ndof);
                     
                titer=toc;
                titertot=titertot+titer;
        
                if printStatus(deltaf, uu_bc(1:Body1.ndof + Body2.ndof), tol, ii, jj, imax, steps, titertot, Gap)
                    break;  
                end  
            end

            if approach == 7 

                lambda_converged = ( norm(lambda_next - lambda) <= tol );             
                lambda = lambda_next;     

            else
                lambda_converged = true;                
            end

        end

    Body1 = SaveResults(Body1,ii,"last"); % options: "all", "last", each by (number) 
    Body2 = SaveResults(Body2,ii,"last");

end
% %##################### Post-Processing ######################
typeM = approachName(approach);
fig_number = 1; 
ShowNodeNumbers = true;
fprintf('Static test, contact approach = %s  \n', typeM);
PrintResults(Body1)
PrintResults(Body2)
gam = 1/pn;
hold on 
vis = "uy"; % options: "frame", "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy" 
Visualization(Body1,fig_number,vis,ShowNodeNumbers);
Visualization(Body2,fig_number,vis,ShowNodeNumbers);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization of contact points and contact method
[ContactPoints, ~] = ContactPointfunc(Body1);  
h1 = plot(ContactPoints(:,1),ContactPoints(:,2),'ok','MarkerFaceColor', 'k', 'MarkerSize', 4);
legend('contact points')
legend(h1, 'contact points'); 
gapStr = sprintf('%.5f', Gap);
fullstr = "Method = " + typeM + ", Total Gap = " + gapStr;
title(fullstr, 'Interpreter', 'latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('total steps is %d  \n', total_steps );
fprintf('total gap is %10.22f  \n', Gap )
