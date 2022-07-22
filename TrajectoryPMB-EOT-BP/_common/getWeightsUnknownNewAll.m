function [weights,updatedExistence] = getWeightsUnknownNewAll(logWeights,oldExistence,skipIndex)

logSumWeights = sum(logWeights,2);

if(skipIndex)
    logSumWeights = logSumWeights - logWeights(:,skipIndex);
end

[log_w,log_sum_w] = normalizeLogWeights(logSumWeights);

aliveProbability = oldExistence*exp(log_sum_w);
updatedExistence = aliveProbability/(aliveProbability+1-oldExistence);
weights = exp(log_w);

end