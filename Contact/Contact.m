function [DofsFunction, Stiffness] = Contact(Body1,Body2,approach)


if approach.Type ~= "None" % we have contact algorithm     
    AimFunction = approach.AimFunction;  
    DofsFunction = AimFunction(Body1,Body2);
    switch approach.perturbation 
        case "automatic"
            DofsFunction_y = @(t,y) AimFunctionWrapper(t,y,Body1,Body2,AimFunction);
            Stiffness = numjac(DofsFunction_y, 0,  [Body1.u(:); Body2.u(:)], DofsFunction, 1e-5, []);
            
        case "incremental"    
            Stiffness = Jac_stepwise(Body1,Body2,AimFunction);
    end       

else
    DOFsNumber = Body1.nx + Body2.nx;
    DofsFunction = zeros(DOFsNumber,1);
    Stiffness = zeros(DOFsNumber);
end


