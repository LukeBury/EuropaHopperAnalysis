function [m2] = R3(m1,w)
m1 = m1';
R = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];
m2 = R*m1;
end