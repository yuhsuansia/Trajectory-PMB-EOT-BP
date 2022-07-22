function [ A, Q ] = getTransitionMatricesNew( scanTime )

A = kron([1 scanTime;0 1],eye(2));
Q = kron([scanTime^3/3 scanTime^2/2;scanTime^2/2 scanTime],eye(2));

end

