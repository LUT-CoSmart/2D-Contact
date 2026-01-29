function Gap = InnerGap(Gapfunc,Body1,Body2)
    
    Gap = 0;
    GapAllPoint = Gapfunc(Body1,Body2);
    for ii = 1:length(GapAllPoint)
        if GapAllPoint(ii) < 0
            Gap = Gap + abs(GapAllPoint(ii));
        end
    end    