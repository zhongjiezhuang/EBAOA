function MTBAOA = jEBAOA(feat,label,opts) %
if isfield(opts,'N'), N = opts.N; end
if isfield(opts,'T'), max_Iter = opts.T; end
if isfield(opts,'rate'), rate = opts.rate; end
if isfield(opts,'kk'), kk = opts.kk; end
if isfield(opts,'k'), k = opts.k; end
if isfield(opts,'fun'), fun = opts.fun; end
dim = size(feat,2);
max_Iter = ceil(dim/50);

fprintf('jEBAOA Working SF %s %s %s %s %s %s',num2str(rate),num2str(kk),num2str(k),opts.p_name,num2str(opts.times_id),num2str(max_Iter));
UB = 1; LB=0; UB_Contiuous = 100; LB_Contiuous = -100;
% max_Iter=110;
opts.no_of_tasks = 4; iter=0; iter_all = N*max_Iter;
% rate=0.5;
%Construction task
pos = 1:dim;
RF     = jffs('rf',feat,label);
PCC    = jffs('pcc',feat,label);
CHI2   = jffs('chi2',feat,label);
NCA    = jffs('nca',feat,label);
FilterSF = [PCC.sf;CHI2.sf;RF.sf;NCA.sf];
FilterS = [PCC.s;CHI2.s;RF.s;NCA.s];
rf_unpos = setdiff(pos,RF.sf(1,1:dim*0.5));     %RF
pcc_unpos = setdiff(pos,PCC.sf(1,1:dim*0.5));   %PCC
chi2_unpos = setdiff(pos,CHI2.sf(1,1:dim*0.5)); %CHI2
nca_unpos = setdiff(pos,NCA.sf(1,1:dim*0.5));   %NCA
Tasks(1).unpos=rf_unpos;     Tasks(1).D=dim - size(rf_unpos,2);      
Tasks(2).unpos=pcc_unpos;    Tasks(2).D=dim - size(pcc_unpos,2);
Tasks(3).unpos=chi2_unpos;   Tasks(3).D=dim - size(chi2_unpos,2);
Tasks(4).unpos=nca_unpos;    Tasks(4).D=dim - size(nca_unpos,2);
%initialization
for i = 1 : floor(N/4)
%     population(i) = Particle();
    population(4*i-3) = initialize(Particle(),kk,PCC.sf,dim,rate);
%     population(i) = initialize3(population(i),dim);
    population(4*i-3).skill_factor=0;
    population(4*i-2) = initialize(Particle(),kk,CHI2.sf,dim,rate);
    population(4*i-2).skill_factor=0;
    population(4*i-1) = initialize(Particle(),kk,RF.sf,dim,rate);
    population(4*i-1).skill_factor=0;
    population(4*i) = initialize(Particle(),kk,NCA.sf,dim,rate);
    population(4*i).skill_factor=0;
end

for i = floor(N/4)*4:N
    population(i) = initialize(Particle(),kk,NCA.sf,dim,rate);
    population(i).skill_factor=0;
end
for i = 1 : N
    population(i) = evaluate(feat,label,population(i),Tasks,opts);
    iter = iter + population(i).calls;
end
factorial_cost=zeros(1,N); no_of_tasks = 4; 
gbest = zeros(1,dim);
gobj = inf;
gerror = inf;
gselefeat = 0;
for i = 1:no_of_tasks      
    for j = 1:N
        factorial_cost(j)=population(j).factorial_costs(i);
    end                     
    [xxx,y]=sort(factorial_cost);
    population=population(y);  % reorder the population according to factorial_cost of current task            
    for j=1:N
        population(j).factorial_ranks(i)=j; 
    end
    Tasks(i).tobj = population(1).factorial_costs(i);
    Tasks(i).tbest = population(1).rnvec;
    Tasks(i).terror =  population(1).factorial_error(i);
    Tasks(i).tselefeat =  population(1).factorial_selefeat(i);
    if Tasks(i).tobj<gobj
        gobj = Tasks(i).tobj;
        gbest(1,:) = Tasks(i).tbest;
        gerror = Tasks(i).terror;
        gselefeat = Tasks(i).tselefeat; 
    end
end
for i=1:N
    [xxx,yyy]=min(population(i).factorial_ranks);
    x=find(population(i).factorial_ranks == xxx);
    equivalent_skills=length(x);
    if equivalent_skills>1  % If having best fitness on multiple tasks, random choose one and set the factorial_costs of others as inf
        population(i).skill_factor=x(1+round((equivalent_skills-1)*rand(1)));
        tmp=population(i).factorial_costs(population(i).skill_factor);
        population(i).factorial_costs(1:no_of_tasks)=inf;
        population(i).factorial_costs(population(i).skill_factor)=tmp;
    else                    % else, just set the skill_factor and set the factorial_costs of others as inf
        population(i).skill_factor=yyy;
        tmp=population(i).factorial_costs(population(i).skill_factor);
        population(i).factorial_costs(1:no_of_tasks)=inf;
        population(i).factorial_costs(population(i).skill_factor)=tmp;
    end
end

MOP_Max=1;MOP_Min=0.2;Alpha=5;Mu=0.499;
curve.cost = zeros(1,max_Iter); curve.error = zeros(1,max_Iter); curve.sele_feat = zeros(1,max_Iter);
curve.cost(1) = gobj; curve.error(1) = gerror; curve.sele_feat(1) = gselefeat;
fprintf('\nIteration %d Best (MTBAOA_SF)= %f',1,curve.cost(1))
t=1;

while t<(max_Iter+1) %&& iter< (iter_all+1)  %Main loop
    MOP=1-((t)^(1/Alpha)/(max_Iter)^(1/Alpha));   % Probability Ratio
    MOA=MOP_Min+t*((MOP_Max-MOP_Min)/max_Iter); %Accelerated function
    
    %Update the Position of solutions
    for i=1:N   % if each of the UB and LB has a just value
        Xnew = population(i);
        skill_factor = population(i).skill_factor;
        for j=1:dim
            r1=rand();
            if r1>MOA
                r2=rand();
                if r2>0.5
                    flag = transfer(((UB_Contiuous-LB_Contiuous)*Mu+LB_Contiuous)/(MOP+eps));
                else
                    flag = transfer(-((UB_Contiuous-LB_Contiuous)*Mu+LB_Contiuous)/(MOP+eps));
                end
                
                if flag ==0
%                     if rand<0.5
                        Xnew.rnvec(1,j)=gbest(1,j)*~population(i).rnvec(1,j);
%                     else
%                         Xnew.rnvec(1,j)=Tasks(skill_factor).tbest(1,j)*~population(i).rnvec(1,j);
%                     end
                else
%                     if rand<0.5
                        Xnew.rnvec(1,j)=gbest(1,j)*population(i).rnvec(1,j);
%                     else
%                         Xnew.rnvec(1,j)=Tasks(skill_factor).tbest(1,j)*population(i).rnvec(1,j);
%                     end
                end
            else
                r3=rand();
                if r3>0.5
                    flag = transfer(-MOP*((UB_Contiuous-LB_Contiuous)*Mu+LB_Contiuous));                  
                else
                    flag = transfer(MOP*((UB_Contiuous-LB_Contiuous)*Mu+LB_Contiuous));
                end
                
                if flag ==0
%                     if rand<0.5
                        Xnew.rnvec(1,j)=gbest(1,j)+~population(i).rnvec(1,j);
%                     else
%                         Xnew.rnvec(1,j)=Tasks(skill_factor).tbest(1,j)+~population(i).rnvec(1,j);
%                     end
                else
%                     if rand<0.5
                        Xnew.rnvec(1,j)=gbest(1,j)+population(i).rnvec(1,j);
%                     else
%                         Xnew.rnvec(1,j)=Tasks(skill_factor).tbest(1,j)+population(i).rnvec(1,j);
%                     end
                end
            end
        end       
        Flag_UB=Xnew.rnvec(1,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew.rnvec(1,:)<LB; % check if they exceed (down) the boundaries
        Xnew.rnvec(1,:)=(Xnew.rnvec(1,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
        
         population(i) = evaluate(feat,label,Xnew,Tasks,opts); iter = iter + 1;
         factorial_cost=zeros(1,N); 
         for ii = 1:no_of_tasks      
            for j = 1:N
                factorial_cost(j)=population(j).factorial_costs(ii);
            end                     
            [xxx,y]=sort(factorial_cost);
            population=population(y);  % reorder the population according to factorial_cost of current task            
            for j=1:N
                population(j).factorial_ranks(ii)=j; 
            end
            Tasks(ii).tobj = population(1).factorial_costs(ii);
            Tasks(ii).tbest = population(1).rnvec;
            Tasks(ii).terror =  population(1).factorial_error(ii);
            Tasks(ii).tselefeat =  population(1).factorial_selefeat(ii);
            if Tasks(ii).tobj<gobj
                gobj = Tasks(ii).tobj;
                gbest(1,:) = Tasks(ii).tbest;
                gerror = Tasks(ii).terror;
                gselefeat = Tasks(ii).tselefeat; 
            end
        end
        for ii=1:N
            [xxx,yyy]=min(population(ii).factorial_ranks);
            x=find(population(ii).factorial_ranks == xxx);
            equivalent_skills=length(x);
            if equivalent_skills>1  % If having best fitness on multiple tasks, random choose one and set the factorial_costs of others as inf
                population(ii).skill_factor=x(1+round((equivalent_skills-1)*rand(1)));
                tmp=population(ii).factorial_costs(population(ii).skill_factor);
                population(ii).factorial_costs(1:no_of_tasks)=inf;
                population(ii).factorial_costs(population(ii).skill_factor)=tmp;
            else                    % else, just set the skill_factor and set the factorial_costs of others as inf
                population(ii).skill_factor=yyy;
                tmp=population(ii).factorial_costs(population(ii).skill_factor);
                population(ii).factorial_costs(1:no_of_tasks)=inf;
                population(ii).factorial_costs(population(ii).skill_factor)=tmp;
            end
        end   
    end
    
    curve.cost(t) = gobj; curve.error(t) = gerror; curve.sele_feat(t) = gselefeat;
    if mod(t,10)==0
        fprintf('\nIteration %d Best (MTBAOA_SF)= %f',t,curve.cost(t))
        toc
    end
%     fnum(t) = length(find(X == 1));
    t = t + 1;
end
% Select features based on selected index
Pos   = 1:dim;
Sf    = Pos(gbest == 1); 
sFeat = feat(:,Sf); 
% Store results
MTBAOA.sf = Sf; 
MTBAOA.ff = sFeat;
MTBAOA.nf = length(Sf);
MTBAOA.c  = curve.cost;
MTBAOA.r  = curve.error;
MTBAOA.selefeat  = curve.sele_feat;
MTBAOA.f  = feat;
MTBAOA.l  = label;
% MTBAOA.fnum = fnum;
MTBAOA.fnum = 1;
MTBAOA.t = t-1;
end

function bin_x=transfer(x)
s=abs(tanh(x));% V2 transfer function
if s > rand()
    bin_x=0;
else
    bin_x=1;
end
end





