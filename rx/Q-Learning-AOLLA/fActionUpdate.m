function [action1,action2] = fActionUpdate(epsilon,Qtable,beta1new,beta2new)
 i = rand(1);
lo4 = find(Qtable(1,:)== beta1new);
lo5 = find(Qtable(1,:)== beta2new);
 if i > epsilon
     action1lo = min(find(Qtable(2:end,lo4) == max(Qtable(2:end,lo4)))+1);
     action1 = Qtable(action1lo,1);
     action2lo = min(find(Qtable(2:end,lo5) == max(Qtable(2:end,lo5)))+1);
     action2 = Qtable(action2lo,1);
 else
     action1lo = randi(3,1)+1;
     action1 = Qtable(action1lo,1);
     action2lo = randi(3,1)+1;
     action2 = Qtable(action2lo,1);
 end
end