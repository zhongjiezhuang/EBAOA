clear, clc, close; warning('off','all');
dataset_num = 24;
times_count = 10; 
result_all = zeros(times_count,dataset_num*4); % 1,2,3,4列为第一个dataset的error,num,fitness,time
result_conclude = zeros(dataset_num*3,4); %第一行是第一个函数20次结果的均值，第二行是std。
Problem  = {'Yale','ORL','colon','warpAR10P','warpPIE10P','PCMAC_S', 'lung','lymphoma','RELATHE_S', 'GLIOMA','Breast3', 'TOX_171',...
      'Prostate_GE','leukemia','ALL_AML','Command','Name','Adenocarcinoma', 'nci9', 'arcene','orlraws10P','CLL_SUB_111','SMK_CAN_187','GLI_85'};

%% MAIN LOOP
for i = 1:length(Problem)
    p_name = Problem{i};        
    traindata = load(['mtdata\',p_name]);
    data = getfield(traindata, 'data');
    feat = data(:,1:end-1); 
    label = data(:,end);  
    for times_id = 1:times_count
        %% Basic settings EC algorithms    
        opts.k = 1; % Number of k in K-nearest neighbor
        opts.N  = 50;%20;     % number of solutions
        opts.thres = 0.6;
        opts.T  = 100;    % number of iterations in each phase
        opts.p_name = p_name;  
        opts.times_id = times_id;
        opts.fun = @jFitnessFunction;
        %% Swarm-based  
        tic
        opts.rate = 0.05 ;opts.kk = 10;Temp = jEBAOA(feat,label,opts); toc 
        result_all(times_id,(4*i-3):4*i) = [Temp.c(:,Temp.t) Temp.r(:,Temp.t) Temp.nf toc];
    end 
    result_mean = mean(result_all,1); 
    result_min = min(result_all,[],1);
    result_std = std(result_all,1);  
    result_conclude(3*i-2,1:4) = result_mean(1,(4*i-3):4*i); %error,number of selected features , fitness value
    result_conclude(3*i-1,1:4) = result_min(1,(4*i-3):4*i);
    result_conclude(3*i,1:4) = result_std(1,(4*i-3):4*i);
end 

