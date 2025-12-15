function Gap = GapNodes(ContactBody,TargetBody)

Gap = 0;
DofsAtNode_cont = ContactBody.DofsAtNode;
ContactNode_cont = ContactBody.contact.nodalid;

% current position of the "possible contact" nodes of the Contact Body
ContactPoints_X = ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,1)) + ...
                  ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,1));   ... % coords on X axis;

ContactPoints_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                   ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));

ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact body of the contact surfaces

%% TODO: make it over all points simultaneously, working with array or in parallel  
for ii = 1:size(ContactPoints,1) % loop over all contact points
  
    ContactPoint = ContactPoints(ii,:);
 
    % Search the attributes of the corresponding point on the target surface 
    Outcome = FindPoint(TargetBody,ContactPoint);

    % Checking the condition of the penalty approach
    if Outcome.Gap < 0 % we have meaningful outcome from the search
        Gap= Gap + abs(Outcome.Gap);   
    end   
end


