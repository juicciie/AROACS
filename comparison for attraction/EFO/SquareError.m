function E = SquareError(C)
% Optimal point is centered to 40 for each parameter
E=sum((C-40).^2);
end
