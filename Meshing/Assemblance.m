function [uu_bc, deltaf,lambda_next] = Assemblance(Body1,Body2,Fc,Kc,GapNab,approach,penalty,lambda)
    
    lambda_next = lambda;
    bc = [Body1.bc Body2.bc]; % total logical vector of constrains
    Fext = [Body1.Fext.vec; Body2.Fext.vec]; % Assemblance of external forces
    
    % Assemblance of elastic stiffness matrices
    Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
          zeros(Body2.nx,Body1.nx)            Body2.Fint.K];
    
    Fe = [Body1.Fint.vec; Body2.Fint.vec]; % Assemblance of elastic forces
    
    ff =  Fe - Fext + Fc;

    % Remuval of the fixed dofs
    Ke_bc = Ke(bc,bc);  
    Kc_bc = Kc(bc,bc);  
    K_bc = Ke_bc + Kc_bc; % stiffness matrices assemblance
    ff_bc = ff(bc);
    GapNab_bc = GapNab(bc);

    if (approach == 5) || (approach == 8) % Lagrange multiplier & Lagrange multiplier (nonlinear constrain)
                                                                                 

        K_bc = Ke_bc + lambda * Kc_bc; % for Lagrange multiplier Kc = zero-matrix

        K_bc = [     K_bc GapNab_bc;
                GapNab_bc'        0];  
        
        ff_bc = [ff_bc + lambda*GapNab_bc; 0];

    elseif (approach == 9) || (approach == 10) % perturbed Lagrangian method & perturbed Lagrangian method (nonlinear constrain)

        K_bc = Ke_bc + lambda * Kc_bc; % for perturbed Lagrangian method, Kc = zero-matrix

        K_bc = [     K_bc GapNab_bc;
                GapNab_bc'     -1/penalty];   
        ff_bc = [ff_bc + lambda*GapNab_bc; lambda/penalty];
            
    elseif approach == 6 % very simplified penalty

        Kc_bc = penalty * (GapNab_bc * GapNab_bc');
        K_bc  = Ke_bc + Kc_bc;
    
    elseif approach == 7 % Augumented Lagrange
                
        ff_bc = ff_bc + lambda; % contributions from the updated contact forces after the previous iteration  
        
    end
   
    uu_bc = -K_bc\ff_bc; 

    if (approach == 5) || (approach == 8) || (approach == 9) || (approach == 10)
       lambda_next = lambda + uu_bc(end);
    elseif approach == 7
        lambda_next = lambda + Fc([Body1.bc Body2.bc]); % update 
    elseif approach == 11
       % if norm(lambda) > 1e-5 
       %  lambda = lambda + uu_bc(end-m+1:end); 
       % end 
    end   

        
    deltaf = ff_bc(1:size(Ke_bc,1))/norm(Fext(bc));