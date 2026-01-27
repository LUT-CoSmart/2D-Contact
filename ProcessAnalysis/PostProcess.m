
fig_number = 1; 
ShowNodeNumbers = true;
fprintf('Static test, contact approach = %s  \n', approach.Name);
PrintResults(Body1)
PrintResults(Body2)

if approach.Type == "Penalty"
    gam = 1/approach.penalty;
else
    gam = NaN;
end

hold on 
vis = "uy"; % options: "frame", "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy" 
Visualization(Body1,fig_number,vis,ShowNodeNumbers);
Visualization(Body2,fig_number,vis,ShowNodeNumbers);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization of contact points and contact method
[ContactPoints, ~] = ContactPointfunc(Body1);  
h1 = plot(ContactPoints(:,1),ContactPoints(:,2),'ok','MarkerFaceColor', 'k', 'MarkerSize', 4);
legend('contact points')
legend(h1, 'contact points'); 
gapStr = sprintf('%.5f', Gap);
fullstr = "Method = " + approach.Name + ", Total Gap = " + gapStr;
title(fullstr, 'Interpreter', 'latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('total steps is %d  \n', total_steps );
fprintf('total gap is %10.22f  \n', Gap )
Gapfunc(Body1,Body2)