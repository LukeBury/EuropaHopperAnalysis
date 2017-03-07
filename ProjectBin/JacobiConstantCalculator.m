% ------------------------------------------------------------------------
%%% Jacobi Constant Calculator
% ------------------------------------------------------------------------
%%% Inputs:
% 1) u (secondary body mass / primary body mass) .[1x1]
% 2) Positions of hopper in BCR frame ............[nx3]
% 3) Velocities of hopper in BCR frame ...........[nx3]
% 4) Position of primary body in BCR frame .......[1x3]
% 5) Position of secondary body in BCR frame .....[1x3]
%%% Outputs:
% 1) Jacobi constants ............................[nx1]
function [JCs] = JacobiConstantCalculator(u,rH,vH,rB1,rB2)
JCs = zeros(length(rH),1);
dimTest = size(rH,2) + size(vH,2) + size(rB1,2) + size(rB2,2);
if dimTest == 12
    for k = 1:length(rH)
        r1 = norm(rB1 - rH(k,:));
        r2 = norm(rB2 - rH(k,:));
        JCs(k) = rH(k,1)^2 + rH(k,2)^2 + 2*(1-u)/r1 + 2*u/r2 - vH(k,1)^2 - vH(k,2)^2 - vH(k,3)^2;
    end
else
    warning('Check input dimensions\n')
end
end