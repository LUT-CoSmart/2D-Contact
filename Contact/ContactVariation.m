function [Fc,Kc,GapNab,Gap] = ContactVariation(Body1,Body2,penalty,approach,ContactPointfunc,Gapfunc)
    
    sqrtEps = sqrt(eps);
    h = 2*sqrtEps;
     
    % Current (!!total!!) gap calculation 
    Gap = Gapfunc(Body1,Body2); 
    
    % Current gap calculation with distribution over all DOFs
    % GapDOFs = GapPointDistribution(Body1,Body2,ContactPointfunc);

    % Backup original coordinates
    u1_backup = Body1.u;
    u2_backup = Body2.u;

    nx = Body1.nx + Body2.nx;
    
    % Initiation 
    Kc = zeros(nx,nx);    
    GapNab = zeros(nx,1);
    I_vec=zeros(nx,1);    
       
    if (approach ~= 5) && (approach ~= 8) && (approach ~= 9) && (approach ~= 10) % all, but Lagrange multiplier
        Fc= ContactForce(Body1,Body2,penalty,approach,ContactPointfunc); % Body1 forces from the projection of Body2          
    else % Lagrange multiplier  methods  
        Fc = zeros(nx,1); % we aren't interested in contact forces
    end
        
    for ii = 1:nx
        % start variationing 
        I_vec(ii)=1;
        
        % this split is to distribute coord. between bodies 
        Body1.u = u1_backup - h*I_vec(1:Body1.nx); 
        Body2.u = u2_backup - h*I_vec(1+Body1.nx:end);

        if (approach ~= 5) &&  (approach ~= 6) && (approach ~= 8) && (approach ~= 9) && (approach ~= 10) % Penalty-, Augmented Lagrange & Nitsche-based approaches

            Fch = ContactForce(Body1,Body2,penalty,approach,ContactPointfunc); % force due to variation
            Kc(:,ii) = (Fc - Fch) / h; 
        
        elseif (approach == 6) || (approach == 5) || (approach == 8) || ... % very simplified Penalty & Lagrange multiplier & Lagrange multiplier (nonlinear constrain)
               (approach == 9) || (approach == 10) % perturbed Lagrangian method & perturbed Lagrangian method (nonlinear constrain)

            Gaph = Gapfunc(Body1,Body2);  % -i          
            GapNab(ii) = (Gap - Gaph)/h;         
                                        
            if (approach == 8) || (approach == 10) 
                % Building Hessian
                I_vec1 = zeros(nx,1); 
                for jj = 1:nx
                    I_vec1(jj) = 1;

                    % this split is to distribute coord. between bodies 
                    Body1.u = u1_backup - h*I_vec(1:Body1.nx) - h*I_vec1(1:Body1.nx); 
                    Body2.u = u2_backup - h*I_vec(1+Body1.nx:end) - h*I_vec1(1+Body1.nx:end);
                    Gaph_mm = Gapfunc(Body1,Body2);  % --   

                    Body1.u = u1_backup - h*I_vec(1:Body1.nx) + h*I_vec1(1:Body1.nx); 
                    Body2.u = u2_backup - h*I_vec(1+Body1.nx:end) + h*I_vec1(1+Body1.nx:end);
                    Gaph_mp = Gapfunc(Body1,Body2);  % -+

                    Body1.u = u1_backup + h*I_vec(1:Body1.nx) - h*I_vec1(1:Body1.nx); 
                    Body2.u = u2_backup + h*I_vec(1+Body1.nx:end) - h*I_vec1(1+Body1.nx:end);
                    Gaph_pm = Gapfunc(Body1,Body2);  % +-

                    Body1.u = u1_backup + h*I_vec(1:Body1.nx) + h*I_vec1(1:Body1.nx); 
                    Body2.u = u2_backup + h*I_vec(1+Body1.nx:end) + h*I_vec1(1+Body1.nx:end);
                    Gaph_pp = Gapfunc(Body1,Body2);  % ++

                    Kc(ii,jj) = 1/(4*h^2) * ( Gaph_mm - Gaph_mp - Gaph_pm + Gaph_pp );

                    I_vec1(ii)=0;                
                end
            end
             

        else   % testing
         
        %     GapDOFsh = GapPointDistribution(Body1,Body2,ContactPointfunc);
        %     Kc(:,ii) = (GapDOFs - GapDOFsh) / h;

        end
  
        % end variationing 
        I_vec(ii)=0;      
    end    

    % returning all back to normal
    Body1.u = u1_backup;
    Body2.u = u2_backup;