function [outputDA] = dataAssociationBP(inputDA)

% perform DA
inputDA = inputDA(2,:)./inputDA(1,:);
sumInputDA = 1 + sum(inputDA,2);
outputDA = 1 ./ (repmat(sumInputDA,[1,size(inputDA,2)]) - inputDA);


% make hard DA decision in case outputDA involves NANs
if(any(isnan(outputDA)))
    outputDA = zeros(size(outputDA));
    [~,index] = max(inputDA(1,:));
    outputDA(index(1)) = 1;
end


end