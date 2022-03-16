function model = SAcosModel()
  model.name = 'SAcos Model';
	model.paramNames = {'kappa_m', 'capacity', 'kappa_r', 'kappa_f', 's'};
	model.lowerbound = [0 1 0 0 0];     % Lower bounds for the parameters
	model.upperbound = [Inf Inf Inf Inf Inf Inf]; % Upper bounds for the parameters
	model.movestd = [0.1, 0.1, 0.1, 0.1, 0.1];
	model.pdf = @SAcosPDF;
	model.start = [5, 2, 5, 5, .2;    % Jm, capacity, kappa_r, Jf, s
                 10, 3, 10, 10, .5;
                 100, 4, 100, 100, 1];   
end