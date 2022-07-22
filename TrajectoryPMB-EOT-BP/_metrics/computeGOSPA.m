function GOSPA = computeGOSPA(targetTracks,targetExtents,estimatedTracks,estimatedExtents,parameters)

numSteps = size(targetTracks,2);
GOSPA.total = zeros(numSteps,1);
GOSPA.loc = zeros(numSteps,1);
GOSPA.miss = zeros(numSteps,1);
GOSPA.false = zeros(numSteps,1);

for n = 1:numSteps
    %extract ground truth
    x_mat.x = zeros(2,0);
    x_mat.e = zeros(2,2,0);
    for i = 1:size(targetTracks,3)
        if ~isnan(targetTracks(1,n,i))
            x_mat.x = [x_mat.x targetTracks(1:2,n,i)];
            x_mat.e = cat(3,x_mat.e,targetExtents(:,:,n,i));
        end
    end
    %extract estimates
    y_mat.x = zeros(2,0);
    y_mat.e = zeros(2,2,0);
    for i = 1:size(estimatedTracks,3)
        if ~isnan(estimatedTracks(1,n,i))
            y_mat.x = [y_mat.x estimatedTracks(1:2,n,i)];
            y_mat.e = cat(3,y_mat.e,estimatedExtents(:,:,n,i));
        end
    end
    [d_gospa,~,decomposed_cost] =  GOSPA_extended(x_mat,y_mat,parameters.GOSPA_p,parameters.GOSPA_c,2);
    GOSPA.total(n) = d_gospa;
    GOSPA.loc(n) = decomposed_cost.localisation;
    GOSPA.miss(n) = decomposed_cost.missed;
    GOSPA.false(n) = decomposed_cost.false;
end

end