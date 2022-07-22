function [trajectoryMetric,X,Y] = computeTrajectoryMetric(targetTracks,targetExtents,estimatedTracks,estimatedExtents,parameters)

X.xState = targetTracks(1:2,:,:);
X.eState = targetExtents;

Y.xState = estimatedTracks(1:2,:,:);
Y.eState = estimatedExtents;

nx = size(targetTracks,3);
X.tVec = zeros(nx,1);
X.iVec = zeros(nx,1);

for i = 1:nx
    X.tVec(i) = find(~isnan(targetTracks(1,:,i)),1,'first');
    X.iVec(i) = find(~isnan(targetTracks(1,:,i)),1,'last')-X.tVec(i)+1;
end

ny = size(estimatedTracks,3);
Y.tVec = zeros(ny,1);
Y.iVec = zeros(ny,1);

for i = 1:ny
    Y.tVec(i) = find(~isnan(estimatedTracks(1,:,i)),1,'first');
    Y.iVec(i) = find(~isnan(estimatedTracks(1,:,i)),1,'last')-Y.tVec(i)+1;
end

[d_traMetric,loc_cost,miss_cost,fa_cost,switch_cost]=LPTrajMetric_cluster_EOT(X,Y,parameters.traMetric_c,parameters.traMetric_p,parameters.traMetric_gamma);
trajectoryMetric.total = d_traMetric;
trajectoryMetric.loc = loc_cost;
trajectoryMetric.miss = miss_cost;
trajectoryMetric.false = fa_cost;
trajectoryMetric.switch = [switch_cost;0];

end