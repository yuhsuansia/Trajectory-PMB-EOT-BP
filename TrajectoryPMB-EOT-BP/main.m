% Yuxuan Xia, 2021
% Implementation of the trajectory PMB filter for extended object tracking
% using belief propagation for the set of all trajectories, including both
% alive and dead trajectories, at each time step. This implementation is
% adapted from the implementation in the following publication:

% F. Meyer and J. L. Williams, “Scalable detection and tracking of geometric extended objects,” 
% IEEE Trans. Signal Process., vol. 69, pp. 6283–6298, Oct. 2021.
% available at https://github.com/meyer-ucsd/EOT-TSP-21

clear variables; close all; clc; addpath ('./_common','./_metrics'); dbstop if error

rng default

% parameters of simulated scenario
numSteps = 100;
numTargets = 10;
meanTargetDimension = 9;
startRadius = 75;
startVelocity = 10;

parameters.numSteps = numSteps;

% main parameters of the statistical model
parameters.scanTime = .2;
parameters.accelerationDeviation = 1;
parameters.survivalProbability = 0.99;
parameters.meanBirths = .01;
parameters.surveillanceRegion = [[-150; 150] [-150; 150]];
parameters.measurementVariance = 1^2;
parameters.meanMeasurements = 5;
parameters.meanClutter = 10;

% prior distribution parameters
parameters.priorVelocityCovariance = diag([10^2;10^2]);
parameters.priorExtent2 = 1000;                       
parameters.priorExtent1 = [[meanTargetDimension 0];[0 meanTargetDimension]]*(parameters.priorExtent2-3);
parameters.degreeFreedomPrediction = 1000;

% sampling parameters
parameters.numParticles = 1000;
parameters.numNewBornParticles = 2000;
parameters.regularizationDeviation = 0;

% detection and pruning parameters
parameters.detectionThreshold = .5;
parameters.thresholdPruning = 10^(-3);
parameters.minimumTrackLength = 1;
parameters.thresholdTruncate = 10^(-4);

% Poisson point process for undetected objects parameters
parameters.PoissonPointProcessEnable = true;

% censoring and measurement reordering parameters
parameters.freeThreshold = 0.9;
parameters.clusterThreshold = 0.9;
parameters.minClusterElements = 1;

% message passing parameters
parameters.numOuterIterations = 3;

% backward simulation parameters
parameters.backwardSimulation = 20;

% GOSPA parameters
parameters.GOSPA_p = 1;
parameters.GOSPA_c = 20;

% LP trajectory metric parameters
parameters.traMetric_p = parameters.GOSPA_p;
parameters.traMetric_c = parameters.GOSPA_c;
parameters.traMetric_gamma = 2;

% generate true start states
[startStates,startMatrixes] = getStartStates(numTargets,startRadius,startVelocity,parameters);
appearanceFromTo = [[3;83],[3;83],[6;86],[6;86],[9;89],[9;89],[12;92],[12;92],[15;95],[15;95]];

% generate true track
[targetTracks,targetExtents] = generateTracksUnknown(parameters,startStates,startMatrixes,appearanceFromTo,numSteps);

% generate measurements
measurements = generateClutteredMeasurements(targetTracks, targetExtents, parameters);

% recursive estimation of the set of all trajectories, TPMB-BP
% [estimatedTracks,estimatedExtents,estimatedParticles,trajectoryMetricFilterAccumulated,runTime] = eotEllipticalShapeRFSofAllTrajectories(measurements,targetTracks,targetExtents,parameters);

% recursive estimation of the set of targets, with trajectory building using meta-data, PMB-BP
[estimatedTracks,estimatedExtents,estimatedParticles,trajectoryMetricFilterAccumulated,runTime] = eotEllipticalShapeRFS(measurements,targetTracks,targetExtents,parameters);

% compute GOSPA for filtering
GOSPA_Filter = computeGOSPA(targetTracks,targetExtents,estimatedTracks,estimatedExtents,parameters);

% compute LP trajectory metric for filtering
[trajectoryMetricFilter,X,Y] = computeTrajectoryMetric(targetTracks,targetExtents,estimatedTracks,estimatedExtents,parameters);

% compute smoothed estimates using backward simulation
[smoothedTracks,smoothedExtents] = computeSmoothedEstimates(estimatedParticles,parameters);

% compute GOSPA for smoothing
GOSPA_Smoother = computeGOSPA(targetTracks,targetExtents,smoothedTracks,smoothedExtents,parameters);

% compute LP trajectory metric for smoothing
[trajectoryMetricSmoother,~,Z] = computeTrajectoryMetric(targetTracks,targetExtents,smoothedTracks,smoothedExtents,parameters);

% show filtering results
showResults(X,Y);title('Filtering result','Interpreter','latex')

% show smoothing results
showResults(X,Z);title('Smoothing result','Interpreter','latex')

