function [averagevelocity] = velocitycal(str1)
averagevelocity = zeros(100,1);

for k =100000:500:150000
    
    data = importfile(strcat(str1, int2str(k), '.liggghts'));
    data = data{:,:};
    

l = length(data);
v = zeros(l,1);
for i =1:l-1
    for j = (i+1):l     
        v(i) = ((data(i,10))^2 + (data(i,11))^2 + (data(i,12))^2)^0.5;   
    end 
end
   averagevelocity(1+(k-100000)/500) = mean(v); 
end
   
end