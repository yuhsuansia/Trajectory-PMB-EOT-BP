% Florian Meyer, 2020

function [matrixesOut] = getSquare2Fast(matrixesIn)

matrixesOut = matrixesIn;

matrixesOut(1,1,:) =  matrixesIn(1,1,:).*matrixesIn(1,1,:) + matrixesIn(2,1,:).*matrixesIn(2,1,:);
matrixesOut(2,1,:) =  matrixesIn(2,1,:).*matrixesIn(1,1,:) + matrixesIn(2,2,:).*matrixesIn(2,1,:);

matrixesOut(1,2,:) =  matrixesIn(1,1,:).*matrixesIn(1,2,:) + matrixesIn(1,2,:).*matrixesIn(2,2,:);
matrixesOut(2,2,:) =  matrixesIn(2,1,:).*matrixesIn(1,2,:) + matrixesIn(2,2,:).*matrixesIn(2,2,:);
   

end

