classdef Particle    
    properties
        rnvec; % (genotype)--> decode to find design variables --> (phenotype)
        factorial_costs;
        factorial_error;
        factorial_selefeat;
        factorial_ranks;
        scalar_fitness;
        skill_factor;
        calls = 0;
        pbestrnvec;
        pbestCosts;
        pbestError;
        pbestSeleFeat;
    end    
    methods        
          function object = initialize(object,k,task,dim,rate)  %轮盘赌
            object.rnvec = zeros(1,dim);     
            object.rnvec(1,RouletteSelection(k,ceil(rate*dim),task)) = 1;           
            object.pbestrnvec = object.rnvec;
          end
       
        function [object] = evaluate(feat,label,object,Tasks,opts) 
            fun = opts.fun;
            no_of_tasks = opts.no_of_tasks;
            if object.skill_factor == 0
                for i = 1:no_of_tasks
                    x = object.rnvec;     x(1,Tasks(i).unpos)=0;       
                    [object.factorial_costs(i),object.factorial_error(i),object.factorial_selefeat(i)] = fun(feat,label,x,opts);  
                    object.pbestCosts(i) = object.factorial_costs(i);
                    object.pbestError(i) = object.factorial_error(i);
                    object.pbestSeleFeat(i) = object.factorial_selefeat(i);
                    object.calls = object.calls + 1;
                end
            else
                object.factorial_costs(1:no_of_tasks)=inf;
                for i = 1:no_of_tasks
                    if object.skill_factor == i
                        x = object.rnvec;     x(1,Tasks(i).unpos)=0;       
                        [object.factorial_costs(i),object.factorial_error(i),object.factorial_selefeat(i)] = fun(feat,label,x,opts);                        
                        object.calls = object.calls + 1;
                        if object.factorial_costs(i)<object.pbestCosts(i)
                           object.pbestCosts(i) = object.factorial_costs(i);
                           object.pbestError(i) = object.factorial_error(i);
                           object.pbestSeleFeat(i) = object.factorial_selefeat(i); 
                        end
                        if object.factorial_costs(i)<Tasks(object.skill_factor).tobj
                            Tasks(object.skill_factor).tobj = object.factorial_costs(i);
                            Tasks(object.skill_factor).tbest = object.rnvec;
                            Tasks(object.skill_factor).terror = object.factorial_error(i);
                            Tasks(object.skill_factor).tselefeat = object.factorial_selefeat(i); 
                        end
                        break;
                    end
                end
            end
        end              
        
        % position update
        function object=positionUpdate(object)                       
%             object.rnvec=object.rnvec+object.velocity;
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;           
        end
        
        % pbest update
        function object=pbestUpdate(object)
              if object.factorial_costs(object.skill_factor)<object.pbestFitness
                  object.pbestFitness=object.factorial_costs(object.skill_factor);
                  object.pbest=object.rnvec;
              end                   
        end 
 
        
    end
end