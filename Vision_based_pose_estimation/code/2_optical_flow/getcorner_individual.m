function [p0,p1,p2,p3,p4] = getcorner_individual(id)
side = 0.152;
gap = 0.152;
gap2 = 0.178;

%need to define p4
% x = i*(side + gap) where i = remainder - tag id when divided by 12
% y = j * (side + gap) where j = quotient - tage id when divided by 12

i = rem(id,12);
j = floor((id)./12);

p4 = [(side+gap)*i;(side + gap)*j];

    if j >= 3 && j < 6
        p4 = p4 + [0; gap2-gap];
    elseif j >= 6 
        p4 = p4 + [0; (gap2-gap)*2];
    end

p1 = p4 + [side;0];
p2 = p4 + [side; side];
p3 = p4 + [0; side];
p0 = p4 + [side/2; side/2];

end