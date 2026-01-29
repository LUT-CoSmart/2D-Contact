function PostProcess(Body1, Body2, ShowVisualization, vis, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc)

fprintf('Static test, contact approach = %s  \n', approach.Name);
PrintResults(Body1)
PrintResults(Body2)

if approach.Type == "Penalty"
    gam = 1/approach.penalty;
else
    gam = NaN;
end

if ShowVisualization
    hold on 
    Visualization(Body1,vis,ShowNodeNumbers);
    Visualization(Body2,vis,ShowNodeNumbers);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization of contact points and contact method
    [ContactPoints, ~] = ContactPointfunc(Body1);  
    h1 = plot(ContactPoints(:,1),ContactPoints(:,2),'ok','MarkerFaceColor', 'k', 'MarkerSize', 4);
    legend('contact points')
    legend(h1, 'contact points');
    Gap =  Gapfunc(Body1,Body2);
    gapStr = sprintf('%.5f', Gap);
    fullstr = "Method = " + approach.Name + ", Total Gap = " + gapStr;
    title(fullstr, 'Interpreter', 'latex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    
