function Chi2 = jChi2(feat,label,opts)
if isfield(opts,'Nf'), num_feat = opts.Nf; end

% Perform Chi2 Test
[idx,scores] = fscchi2(feat,label);
% model    = fscnca(feat,label,'Solver','sgd');

% Select features based on selected index
Sf            = idx(1:num_feat);
sFeat         = feat(:, Sf); 
% Store results
Chi2.sf = Sf; 
Chi2.ff = sFeat;
Chi2.nf = num_feat; 
Chi2.f  = feat;
Chi2.l  = label;
Chi2.s  = scores;
end

