% Filter Feature Selection Toolbox 

function model = jffs(type,feat,label)
opts.K  = 3; opts.Nf = size(feat,2);
switch type
  case'nca'    ; fun = @jNeighborhoodComponentAnalysis;
  case'pcc'    ; fun = @jPearsonCorrelationCoefficient;
  case'tv'     ; fun = @jTermVarianceFeatureSelection; 
  case'rf'     ; fun = @jReliefF; 
  case'chi2'   ; fun = @jChi2; 
end
tic;
model = fun(feat,label,opts); 
% Computational time
t = toc;

model.t = t;
% fprintf('\n Processing Time (s): %f % \n',t); fprintf('\n');
end


