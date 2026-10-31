function [strategySeq] = Roulette(P)
%轮盘赌算法，随机num次根据P概率给出最后的策略序号
            P = 1./P;
            [M,N] = size(P);
            P = P ./ sum(P,2);
            P = P * triu(ones(N)); %运用上三角函数将计算P的累计概率
            strategySeq=zeros(M,1);
            for i = 1:M
                strategySeq(i,:) = find(P(i,:) >= rand(),1);
            end
%             [strategySeq,strategySeq2] = find(P >= rand(M,N),1);
end

