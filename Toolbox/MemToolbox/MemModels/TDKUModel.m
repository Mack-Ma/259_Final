
function model = TDKUModel()
 
%   This model contains Target, non-target, distractor, and uniform
%   guessing.
%   The precision of target and non-target items are the same.   K (kappa)
%   The precision of target and distracting items are not the same.   KD
%   
%   Created by Yixuan Ku June 2016  added distraction with a different precision (kappa)

    model.name = 'Target, plus Distractor with different Kappa, plus Uniform Model';

    model.paramNames = {'Pn', 'Pd', 'Pu', 'K', 'Kd'};

    model.lowerbound = [0 0 0 0 0];
	model.upperbound = [1 1 1 Inf Inf];

    model.pdf = @TDKUModelPDF;
    model.generator = @TDKUModelGenerator;
	model.start = [0.1, 0.1, 0.2, 32, 33;  % 'Pn', 'Pd', 'Pu', 'K', 'Kd'
                 0.1, 0.1, 0.3, 15, 16;
                 0.3, 0.3, 0.1, 8, 9];
%     model.start = [0 0 0.209 6.356 7465735704955300;
%         0 0 0.209 6.356 7465735704955311]; %% sometimes Kd is a huge number
             
	model.movestd = [0.02, 0.02, 0.02 0.1 0.1];
    model.modelPlot = @nontarget_plot;
%     model.modelPlot = @distractor_plot;
    
    % Example of a possible .priorForMC:
  %   model.priorForMC = @(p) (betapdf(p(1),1.25,2.5) * ... % for Pu
  %   betapdf(p(2),1.25,2.5) * ... % for Pn
  %   lognpdf(deg2k(p(3)),2,0.5)); % for sd

end

% Use customed Plot to make a plot of errors centered on non-target distractors 
function figHand = nontarget_plot(data, params, varargin)
    d.errors = [];
    for i=1:length(data.errors)
      d.errors = [d.errors; circdist(data.errors(i), data.distractors(:,i))];
    end
    if isstruct(params) && isfield(params, 'vals')
      params = MCMCSummarize(params, 'maxPosterior');
    end
    m = StandardMixtureModel();
    f = [1-(params(1)/size(data.distractors,1)) params(4)];
    figHand = PlotModelFit(m, f, d, 'NewFigure', true, 'ShowNumbers', false);
    title('Error relative to nontarget locations', 'FontSize', 14);
    topOfY = max(ylim); 
    txt = sprintf('Pn: %.3f\nsd: %0.2f\n', params(1), k2sd(params(4))*180/pi);
    text(180, topOfY-topOfY*0.05, txt, 'HorizontalAlignment', 'right');
end

% Use customed Plot to make a plot of errors centered on ku distractors 
function figHand = distractor_plot(data, params, varargin)
    d.errors = [];
    for i=1:length(data.errors)
      d.errors = [d.errors; circdist(data.errors(i), data.kudistractors(:,i))];
    end
    if isstruct(params) && isfield(params, 'vals')
      params = MCMCSummarize(params, 'maxPosterior');
    end
    m = StandardMixtureModel();
    f = [1-(params(2)/size(data.kudistractors,1)) params(5)];
    figHand = PlotModelFit(m, f, d, 'NewFigure', true, 'ShowNumbers', false);
    title('Error relative to distractor locations', 'FontSize', 14);
    topOfY = max(ylim); 
    txt = sprintf('Pd: %.3f\nsd: %0.2f\n', params(2), k2sd(params(5))*180/pi);
    text(180, topOfY-topOfY*0.05, txt, 'HorizontalAlignment', 'right');
end

% The model's likelihood function.
function p = TDKUModelPDF(data, Pn, Pd, Pu, K, Kd)
  % Parameter bounds check
  if Pn+Pd+Pu > 1
    p = zeros(size(data.errors));
    return;
  end
  
  if(~isfield(data, 'distractors'))
    error('The swap model requires that you specify the non-target.')
  end
  
  if(~isfield(data, 'kudistractors'))
    error('The swap model requires that you specify the ku distractors.')
  end
  
  nDistractors = size(data.distractors,1);
  nKUDistractors = size(data.kudistractors,1);
  p = (1-Pn-Pd-Pu).*vonmisespdf(data.errors(:),0,K) + ...
          (Pu).*unifpdf(data.errors(:), -180, 180);
        
  % Allow for the possibility of NaN's in distractors, as in the case where
  % different trials have different set sizes
  numDistractorsPerTrial = sum(~isnan(data.distractors));
  for i=1:nDistractors
    pdfOut = vonmisespdf(data.errors(:), data.distractors(i,:)', K);
    pdfOut(isnan(pdfOut)) = 0;
    p = p + (Pn./numDistractorsPerTrial(:)).*pdfOut; 
  end
  
  numKUDistractorsPerTrial = sum(~isnan(data.kudistractors));
  for i=1:nKUDistractors
    pdfOut = vonmisespdf(data.errors(:), data.kudistractors(i,:)', Kd);
    pdfOut(isnan(pdfOut)) = 0;
    p = p + (Pn./numKUDistractorsPerTrial(:)).*pdfOut; 
  end
  
end

% TDKU model random number generator
function y = TDKUModelGenerator(params,dims,displayInfo)
  n = prod(dims);

  % Assign types to trials, 0 for remembered, -1 for guess, 1* for
  % non-target (11, 12, ...), 2* for ku distractors (21, 22, ...)  
  r = rand(n,1);
  which = zeros(n,1); % default = remembered
  numDistractorsPerTrial = sum(~isnan(displayInfo.distractors),1)';
  numKUDistractorsPerTrial = sum(~isnan(displayInfo.kudistractors),1)';
  which(r<params{1}+params{2}+params{3}) = 10+ceil(rand(sum(r<params{1}+params{2}+params{3}), 1) ...
    .* numDistractorsPerTrial(r<params{1}+params{2}+params{3})); % swap to random distractor, 11,12,...
  which(r<params{2}+params{3}) = 20+ceil(rand(sum(r<params{2}+params{3}), 1) ...
    .* numKUDistractorsPerTrial(r<params{2}+params{3})); % swap to ku distractor, 21,22,...
  which(r<params{3}) = -1; % guess
  
  % Fill in with errors
  y = zeros(n,1);
  y(which==-1) = rand(sum(which==-1), 1)*360 - 180;   % guessing
  y(which==0)  = vonmisesrnd(0,params{4}, [sum(which==0) 1]);  % remembered
  
  for d=1:size(displayInfo.distractors,1)
    y(which==d+10) = vonmisesrnd(displayInfo.distractors(d,(which==d+10))', ...
      params{4}, [sum(which==d+10) 1]);     % non-target, with kappa K
  end
  
  for dd=1:size(displayInfo.kudistractors,1)
    y(which==dd+20) = vonmisesrnd(displayInfo.kudistractors(dd,(which==dd+20))', ...
      params{5}, [sum(which==dd+20) 1]);    % distractor, with kappa Kd
  end
  
  % Reshape
  y = reshape(y,dims);
end

