function Bernoulli = updateParticlesNewBernoulliAllTrajectory(Bernoulli,logWeights)

logWeights = sum(logWeights,2);
[log_w,log_sum_w] = normalizeLogWeights(logWeights);

aliveUpdate = exp(log_sum_w);
if(isinf(aliveUpdate))
    Bernoulli.existence = 1;
else
    alive = Bernoulli.existence*aliveUpdate;
    dead = (1-Bernoulli.existence);
    Bernoulli.existence = alive/(dead+alive);
end

if(Bernoulli.existence ~= 0)
    Bernoulli.particleWeights = exp(log_w);
else
    Bernoulli.particleWeights = nan(size(logWeights));
    Bernoulli.particleStartTime = nan(size(Bernoulli.particleStartTime));
    Bernoulli.particleEndTime = nan(size(Bernoulli.particleEndTime));
    Bernoulli.particlesKinematic = nan(size(Bernoulli.particlesKinematic));
    Bernoulli.particlesExtent = nan(size(Bernoulli.particlesExtent));
end

end