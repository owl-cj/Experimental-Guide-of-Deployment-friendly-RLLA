function [Qtable,QvalueNew1,QvalueNew2,lo1,lo2,lo11,lo3] = fQtableCal (reward,Qtable,action1,action2,beta1,beta2,beta1new,beta2new,gama,alpha1,alpha2)
lo1 = find(Qtable(:,1)== action1);
lo11 = find(Qtable(:,1)== action2);
lo2 = find(Qtable(1,:)== beta1);
lo3 = find(Qtable(1,:)== beta2);
lo4 = find(Qtable(1,:)== beta1new);
lo5 = find(Qtable(1,:)== beta2new);

if Qtable (lo1,lo2) == 10000
    Qtable (lo1,lo2) = reward/(1-gama);
else
    %Q score update 还是无穷大的值把他变成0计算Q score， 但Qtable里的值不要变
    Qtable1 = Qtable(:,lo4);
    for i = 2:length(Qtable1)
        if Qtable1(i) == 10000
            Qtable1(i) = 0;
        end
    end            
     Qtable(lo1,lo2) = Qtable(lo1,lo2) + alpha1*(reward+gama*max(Qtable1(2:end)) - Qtable(lo1,lo2));
     QvalueNew1 = Qtable(lo1,lo2);
end
if Qtable (lo11,lo3) == 10000
    Qtable (lo11,lo3) = reward/(1-gama);
else
    %Q score update
    Qtable2 = Qtable(:,lo5);
    for i = 2:length(Qtable2)
        if Qtable2(i) == 10000
            Qtable2(i) = 0;
        end
    end
    Qtable(lo11,lo3) = Qtable(lo11,lo3) + alpha2*(reward+gama*max(Qtable2(2:end)) - Qtable(lo11,lo3));
    QvalueNew2 = Qtable(lo11,lo3);
end

end