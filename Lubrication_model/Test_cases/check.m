clear all 
close all
clc


for i = 100:100:50000

    data = importfile_lammps( strcat('./post/Test',int2str(i),'.lammps'));
    data2 = importfile_lammps( strcat('./post3/Test',int2str(i),'.liggghts'));
    p1lammps(i/100,:) = data{1,:};
    p2lammps(i/100,:) = data{2,:};
    p1liggghts(i/100,:) = data2{1,:};
    p2liggghts(i/100,:) = data2{2,:};

end


figure(1)
hold on
title("particle1 x coordinates") 
plot(p1lammps(:,3)) 
plot(p1liggghts(:,3))
legend('lammps', 'liggghts')

figure(2)
hold on 
title("particle1 y coordinates") 
plot(p1lammps(:,4))
hold on 
plot(p1liggghts(:,4))
legend('lammps', 'liggghts')

figure(3)
hold on 
title("particle1 vx ") 
plot(p1lammps(:,9))
hold on 
plot(p1liggghts(:,9))
legend('lammps', 'liggghts')

figure(4)
hold on 
title("particle1 vy ") 
plot(p1lammps(:,10))
hold on 
plot(p1liggghts(:,10))
legend('lammps', 'liggghts')

figure(5)
hold on 
title("particle1 fx ") 
plot(p1lammps(:,12))
hold on 
plot(p1liggghts(:,12))
legend('lammps', 'liggghts')

figure(6)
hold on 
title("particle1 fy ") 
plot(p1lammps(:,13))
hold on 
plot(p1liggghts(:,13))
legend('lammps', 'liggghts')

figure(7)
hold on 
title("particle1 tz") 
plot(p1lammps(:,21))
hold on 
plot(p1liggghts(:,21))
legend('lammps', 'liggghts')

figure(8)
hold on 
title("particle1 omegaz") 
plot(p1lammps(:,17))
hold on 
plot(p1liggghts(:,17))
legend('lammps', 'liggghts')



