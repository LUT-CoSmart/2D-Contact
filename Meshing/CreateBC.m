
function Body = CreateBC(Body)

    DofsAtNode = Body.DofsAtNode;
    
    % Set vector of linear constraints
    bc = true(1,Body.nx);
    
    NodalIDAt = FindGlobNodalID(Body.P0,Body.loc,Body.shift); % let's find global nodal ID's x=0
    
    bcInd=xlocChosen(DofsAtNode,NodalIDAt,1:2);  % all fixed at clambed end
    if bcInd~=0
        bc(bcInd)=0;                       % number of degrees of freedom of system after linear constraints  
    end
    
    Body.bc = bc;
    Body.ndof = sum(bc);     % Number of unconstrained DOFs  