function approach = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc, Gapfunc, GapfuncPairs)

    approach.penalty = 1e10;  
    
    % sanity check, that approach and subtype are correlating 
    if approachBasis == "Lagrange"

       approach.perturbation = "incremental";  

       allowed = ["Lagrange","perturbed Lagrange"];  

       if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Lagrange, substituted to Lagrange'");
          approachSubtype = "Lagrange";
       end

       approach.Name = approachSubtype; 
       % AimFunction = Gapfunc;
       AimFunction = GapfuncPairs; 

    elseif approachBasis == "Penalty"

        approach.perturbation = "automatic";

        allowed = ["Penalty", "Nitshe-linear", "Nitshe-nonlinear", "Nitshe-nonlinear-all", "Augumented Lagrange"];
        
        if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Penalty, substituted to Penalty'");
          approachSubtype = "Penalty";
        end

        if approachSubtype == "Augumented Lagrange"
            approach.penalty = 1e9;  % decreasing parameter for better stability, method operates with any small penalty 
        end

        approach.Name = approachSubtype; 
        AimFunction = @(Body1,Body2) ContactForce(Body1,Body2,approach,ContactPointfunc);
        
     else
        approach.Type = "None";
        warning("Contact is not activated");
        AimFunction = "None";
        approachBasis = "None";    
    end  

    approach.Type = approachBasis;   
    approach.AimFunction = AimFunction;