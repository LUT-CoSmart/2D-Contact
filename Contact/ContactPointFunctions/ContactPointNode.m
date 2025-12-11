function [ContactPoints, ContactPointsElements] = ContactPointNode(ContactBody)
    
    nloc_cont = ContactBody.nloc;
    DofsAtNode_cont = ContactBody.DofsAtNode;
    ContactNode_cont = ContactBody.contact.nodalid;
    
    % current position of the "possible contact" nodes of the Contact Body
    ContactPoints_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                      ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;
    
    ContactPoints_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                       ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));
    
    ContactPoints = [ContactPoints_X ContactPoints_Y]; % nodes of the contact surfaces of the contact body     
    
    
    ContactPointsElements = ContactBody.contact.nodalid(1); % taking the first node number

    for ii = 2:size(ContactPoints,1) % this works only due to the fact that, if you are ath this function then cont point == cont nodes 
        % idea that on the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc_cont == ii-1, 2) & any(nloc_cont == ii, 2)); 
        ContactPointsElements = [ContactPointsElements; ElemenNumber]; % multiplication to correlate with points
    end 
    
    