clear all 
close all
clc

for i = 100:100:50000

    data = importfile_lammps( strcat('./post/Test',int2str(i),'.lammps')); % lammps hooke
    data2 = importfile_lammps( strcat('./post3/Test',int2str(i),'.liggghts')); % liggghts hooke
    
    
    p1swelling_hard(i/100,:) = data{1,:};
    p2swelling_hard(i/100,:) = data{2,:};
    p1hertz_hard(i/100,:) = data2{1,:};
    p2hertz_hard(i/100,:) = data2{2,:};
   
    
    
end


figure(1)
hold on
title("particle1 x coordinates") 
plot(p1swelling_hard(:,3), 'o','MarkerSize',1) 
plot(p1hertz_hard(:,3),'MarkerSize',2)

figure(2)
hold on 
title("particle1 y coordinates") 
plot(p1swelling_hard(:,4),'o','MarkerSize',1) 
plot(p1hertz_hard(:,4),'MarkerSize',2)


figure(3)
hold on 
title("particle1 vx ") 
plot(p1swelling_hard(:,9),'o','MarkerSize',1) 
plot(p1hertz_hard(:,9),'MarkerSize',2)


figure(4)
hold on 
title("particle1 vy ") 
plot(p1swelling_hard(:,10),'o','MarkerSize',1) 
plot(p1hertz_hard(:,10),'MarkerSize',2)

figure(5)
hold on 
title("particle1 fx ") 
plot(p1swelling_hard(:,12),'o','MarkerSize',1) 
plot(p1hertz_hard(:,12),'MarkerSize',2)



figure(6)
hold on 
title("particle1 fy ") 
plot(p1swelling_hard(:,13),'o','MarkerSize',1) 
plot(p1hertz_hard(:,13),'MarkerSize',2)


figure(7)
hold on 
title("particle1 tz") 
plot(p1swelling_hard(:,21),'o','MarkerSize',1) 
plot(p1hertz_hard(:,21),'MarkerSize',2)


figure(8)
hold on 
title("particle1 omegaz") 
plot(p1swelling_hard(:,17),'o','MarkerSize',1) 
plot(p1hertz_hard(:,17),'MarkerSize',2)

