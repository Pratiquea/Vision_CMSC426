%%% RANSAC
%loading given data set into variables
data_1 = load('hw1_data/data1.mat');
data_2 = load('hw1_data/data2.mat');
data_3 = load('hw1_data/data3.mat');
%splitting data_1 into vectors
X_1 = data_1.pts(1,:).';
Y_1 = data_1.pts(2,:).';
%splitting data_2 into vectors
X_2 = data_2.pts(1,:).';
Y_2 = data_2.pts(2,:).';
%splitting data_3 into vectors
X_3 = data_3.pts(1,:).';
Y_3 = data_3.pts(2,:).';

%Selecting random minimum points from datasets



















function [M, inliers] = ransac(data_1, fittingfn, distfn, degenfn, s, t, feedback, ...
 maxDataTrials, maxTrials)

    % Test number of parameters
    error ( nargchk ( 6, 9, nargin ) );
    
    if nargin < 9; maxTrials = 1000;    end;
      if nargin < 8; maxDataTrials = 100; end;
        if nargin < 7; feedback = 0;        end;

          [rows, npts] = size(data_1);
          
    p = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers (probably should be a parameter)

    bestM = NaN;      % Sentinel value allowing detection of solution failure.
    trialcount = 0;
    bestscore =  0;
    N = 1;            % Dummy initialisation for number of trials.
    
    while N > trialcount

        % Select at random s datapoints to form a trial model, M.
        % In selecting these points we have to check that they are not in
        % a degenerate configuration.
        degenerate = 1;
        count = 1;
        while degenerate
            % Generate s random indicies in the range 1..npts
            % (If you do not have the statistics toolbox with randsample(),
            % use the function RANDOMSAMPLE from my webpage)
            if ~exist('randsample', 'file')
              ind = randomsample(npts, s);
            else
              ind = randsample(npts, s);
            end

            % Test that these points are not a degenerate configuration.
            degenerate = feval(degenfn, data_1(:,ind));
            
            if ~degenerate
                % Fit model to this random selection of data points.
                % Note that M may represent a set of models that fit the data in
                % this case M will be a cell array of models
                M = feval(fittingfn, data_1(:,ind));
                
                % Depending on your problem it might be that the only way you
                % can determine whether a data set is degenerate or not is to
                % try to fit a model and see if it succeeds.  If it fails we
                % reset degenerate to true.
                if isempty(M)
                  degenerate = 1;
                end
              end
              
            % Safeguard against being stuck in this loop forever
            count = count + 1;
            if count > maxDataTrials
              warning('Unable to select a nondegenerate data set');
              break
            end
          end
          
        % Once we are out here we should have some kind of model...
        % Evaluate distances between points and model returning the indices
        % of elements in data_1 that are inliers.  Additionally, if M is a cell
        % array of possible models 'distfn' will return the model that has
        % the most inliers.  After this call M will be a non-cell object
        % representing only one model.
        [inliers, M] = feval(distfn, M, data_1, t);
        
        % Find the number of inliers to this model.
        ninliers = length(inliers);
        
        if ninliers > bestscore    % Largest set of inliers so far...
            bestscore = ninliers;  % Record data for this model
            bestinliers = inliers;
            bestM = M;
            
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/npts;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
          end
          
          trialcount = trialcount+1;
          if feedback
            fprintf('trial %d out of %d         \r',trialcount, ceil(N));
          end

        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
          warning( ...
            sprintf('ransac reached the maximum number of %d trials',...
              maxTrials));
          break
        enmd
      end
      
      if feedback, fprintf('\n'); end

    if ~isnan(bestM)   % We got a solution
      M = bestM;
      inliers = bestinliers;
    else
      M = [];
      inliers = [];
      warning('ransac was unable to find a useful solution');
    end
    












