function [epsilon,lambda,E] = plot_lambdaepsilon(N) 

file=strcat('PCOS_lambdaepsilon_',num2str(N));

run(file);

figure

plot(epsilon,lambda)

xlabel('\epsilon')
ylabel('\lambda')

figure

plot(epsilon,E(:,1),epsilon,E(:,2))

xlabel('\epsilon')
ylabel('E')
legend('E^f','E^s');
end