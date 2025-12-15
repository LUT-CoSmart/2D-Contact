function approach = ApproachSettings(approachBasis,approachSubtype)

    % 
    approach.penalty = 1e10;
    approach.lambda.meaning = 0;
    approach.lambda.converge = false;
    
    % sanity check, that approach and subtype are correlating 
    if approachBasis == "Lagrange"
       allowed = ["Lagrange","Lagrange-nonlinear", "perturbed Lagrangian","perturbed Lagrangian-nonlinear"];  

       if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Lagrange, substituted to Lagrange'");
          approachSubtype = "Lagrange";
       end

    elseif approachBasis == "Penalty"
        allowed = ["Penalty", "Nitshe-linear", "Nitshe-nonlinear", "Nitshe-nonlinear-all", "Augumented Lagrange"];
        
        if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Penalty, substituted to Penalty'");
          approachSubtype = "Penalty";
        end

    end  
   
    approach.Type = approachBasis;
    approach.Name = approachSubtype;    