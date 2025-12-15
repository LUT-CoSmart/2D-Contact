function [uu_bc, deltaf, lambda] = Assemblance(Body1,Body2,Fc,Kc,GapNab,approach)
    
    Name = approach.Name;
    lambda = approach.lambda.meaning; % keeping just for standard Lagrange method
    penalty = approach.penalty;

    bc = [Body1.bc Body2.bc]; % total logical vector of constrains
    Fext = [Body1.Fext.vec; Body2.Fext.vec]; % Assemblance of external forces
    
    % Assemblance of elastic stiffness
    Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
          zeros(Body2.nx,Body1.nx)            Body2.Fint.K];
    
    Fe = [Body1.Fint.vec; Body2.Fint.vec]; % Assemblance of elastic forces
    
    ff =  Fe - Fext + Fc; % residuals

    % Removal of the fixed dofs
    Ke_bc = Ke(bc,bc);  
    Kc_bc = Kc(bc,bc);  
    ff_bc = ff(bc);
    GapNab_bc = GapNab(bc);

    % Methods
    switch Name

        case {"Lagrange", "Lagrange-nonlinear"}
            K_bc = Ke_bc + lambda * Kc_bc; % for Lagrange multiplier Kc is zero-matrix
            K_bc = [     K_bc GapNab_bc;
                    GapNab_bc'        0];  
        
            ff_bc = [ff_bc + lambda*GapNab_bc; 0];    

        case {"perturbed Lagrangian", "perturbed Lagrangian-nonlinear"}             
            K_bc = Ke_bc + lambda * Kc_bc; % for perturbed Lagrangian method is Kc = zero-matrix
    
            K_bc = [     K_bc GapNab_bc;
                    GapNab_bc'     -1/penalty];   
            ff_bc = [ff_bc + lambda*GapNab_bc; lambda/penalty];
        
        case "Augumented Lagrange"
            ff_bc = ff_bc + lambda; % contributions from the updated contact forces after the previous iteration
            lambda = lambda + Fc([Body1.bc Body2.bc]); 
            K_bc = Ke_bc + Kc_bc; % standard stiffness matrices assemblance

        otherwise
            K_bc = Ke_bc + Kc_bc; % standard stiffness matrices assemblance

    end    
   
    uu_bc = -K_bc\ff_bc; 
    deltaf = ff_bc(1:size(Ke_bc,1))/norm(Fext(bc));