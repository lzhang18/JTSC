p = 285;
gap_weekend = 96*2;
gap_weekday = 96*5;
gap_week = 96*7;

weekend_idx = [];
for i=1:52
    weekend_idx = [weekend_idx p:p+gap_weekend];
    p = p + gap_week;
end



    
    
