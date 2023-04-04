function reward = fRewardCal (rewardState)
Nbyte = 0;
PacketDuration(:,1) = [1348;708;484;388;276;228;196;196;164;164;148;132];               %不同帧聚合等级,不同MCS下包长统计
PacketDuration(:,2) = [3908;1988;1348;1028;708;548;484;452;388;356;324;292];
PacketDuration(:,3) = [0;3908;2628;1988;1348;1028;916;836;708;644;580;516];
PacketDuration(:,4) = [0;0;3908;2948;1988;1508;1348;1220;1028;932;836;756];
PacketDuration(:,5) = [0;0;5188;3908;2628;1988;1764;1604;1348;1220;1092;980];
 Tprotocol = 84+10;
 for i = 1:length(rewardState(1,:))-1                                                      %延时统计
    if rewardState(3,i) == 12
      Tdata(i) = PacketDuration(1+rewardState(2,i),5);
    end
    if rewardState(3,i) == 9
      Tdata(i) = PacketDuration(1+rewardState(2,i),4); 
    end
    if rewardState(3,i) == 6
      Tdata(i) = PacketDuration(1+rewardState(2,i),3); 
    end    
    if rewardState(3,i) == 3
      Tdata(i) = PacketDuration(1+rewardState(2,i),2); 
    end
    if rewardState(3,i) == 1
      Tdata(i) = PacketDuration(1+rewardState(2,i),1); 
    end    
end
TdataSum = sum(Tdata)+ length(Tdata)*Tprotocol;
for i = 2:length(rewardState(1,:))   
     if rewardState(1,i) == 1
         Nbyte = Nbyte + 2304*rewardState(3,i-1);
     end
end
Rbyte = Nbyte/TdataSum;
reward = Rbyte;
end