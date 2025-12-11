function [ContactPoints, ContactPointsElements] = ContactPointLinSpace(ContactBody,n)

    nloc_cont = ContactBody.nloc;
    DofsAtNode_cont = ContactBody.DofsAtNode;
    ContactNode_cont = ContactBody.contact.nodalid;
    
    % current position of the "possible contact" nodes of the Contact Body
    ContactNodesPosition_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                      ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;
    
    ContactNodesPosition_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                       ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));
    
    ContactNodesPosition = [ContactNodesPosition_X ContactNodesPosition_Y]; % nodes of the contact surfaces of the contact body     
    
    ContactPoints = ContactNodesPosition(1,:); % taking the first node position, later it will be omitted 
    ContactPointsElements = ContactBody.contact.nodalid(1); % taking the first node number
  
    % transition to Linspace points  
    for ii = 2:size(ContactNodesPosition,1)
        % taking two consecutive nodes
        a = ContactNodesPosition(ii-1,:);
        b = ContactNodesPosition(ii,:);    
        % idea that on the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc_cont == ii-1, 2) & any(nloc_cont == ii, 2)); 
        xx = geospace(a(1),b(1),n)'; % split in x- axis
        yy = geospace(a(2),b(2),n)'; % split in y- axis
    
        if abs(a(1) - b(1)) < sqrt(eps)
            xx = a(1)*ones(n,1);
        elseif abs( a(2) - b(2) ) < sqrt(eps)
            yy = a(2)*ones(n,1); 
        end
    
        % Following the Intercept theorem (см. теорему Фалеса):
        % the function devide x- and y- axes in the same propotions 
        ContactPoints = [ContactPoints; xx(2:end) yy(2:end)]; % excluding the node from the previous devision
        ContactPointsElements = [ContactPointsElements; ElemenNumber*ones(n-1,1)]; % multiplication to correlate with points
    end    