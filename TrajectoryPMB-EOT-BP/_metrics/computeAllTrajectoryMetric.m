function [d_traMetric,loc_cost,miss_cost,fa_cost,switch_cost] = computeAllTrajectoryMetric(targetTracks,targetExtents,estimatedTracks,estimatedExtents,step,parameters)

indexValidTrack = find(~all(isnan(targetTracks(1,1:step,:)),2));

X.xState = targetTracks(1:2,1:step,indexValidTrack);
X.eState = targetExtents(:,:,1:step,indexValidTrack);

nx = length(indexValidTrack);
X.tVec = zeros(nx,1);
X.iVec = zeros(nx,1);
for i = 1:nx
    X.tVec(i) = find(~isnan(X.xState(1,:,i)),1,'first');
    X.iVec(i) = find(~isnan(X.xState(1,:,i)),1,'last')-X.tVec(i)+1;
end

indexValidTrack = find(~all(isnan(estimatedTracks(1,1:step,:)),2));

Y.xState = estimatedTracks(1:2,1:step,indexValidTrack);
Y.eState = estimatedExtents(:,:,1:step,indexValidTrack);

ny = length(indexValidTrack);
Y.tVec = zeros(ny,1);
Y.iVec = zeros(ny,1);

for i = 1:ny
    Y.tVec(i) = find(~isnan(Y.xState(1,:,i)),1,'first');
    Y.iVec(i) = find(~isnan(Y.xState(1,:,i)),1,'last')-Y.tVec(i)+1;
end

% rho = parameters.traMetric_rho;
% time_weights1 = (1-rho)/(1-rho^step)*rho.^(step-(1:step))';
% time_weights2 = time_weights1(2:end);

[d_traMetric,loc_cost,miss_cost,fa_cost,switch_cost] = LPTrajMetric_cluster_EOT(X,Y,parameters.traMetric_c,parameters.traMetric_p,parameters.traMetric_gamma);
d_traMetric = d_traMetric/step;
loc_cost = loc_cost/step;
miss_cost = miss_cost/step;
fa_cost = fa_cost/step;
switch_cost = switch_cost/step;

end