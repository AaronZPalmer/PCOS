function [lambdaf, lambdas,Rf,Rs,Phi0,U,S0f,S1f,S0s,S1s,x,t] = plot_example(num)

file=strcat('PCOS_ex',num2str(num));

run(file);

x=(1:399)'/400;

if num==1
    t=(0:100000)/100000;
    T1=90000*ones(399,1);
    T0=9000*ones(399,1);
else
    t=(0:20000)/20000;
    T1=20000*ones(399,1);
    T0=2000*ones(399,1);
end

figure

plot(x,Rf,x,Rs)

xlabel('x');
ylabel('R');

legend('R_f','R_s');

figure

I0=[];
I1=[];

for xi=1:399
    if S0f(xi)<=T0(xi)
        I0=[I0; xi];
    end
    if S1f(xi)<=T1(xi)
        I1=[I1; xi];
    end
end

fill([x(I0);flip(x(I0))],[S0f(I0)./T1(I0); T0(I0)./T1(I0)],'g');
hold on
fill([x(I1); flip(x(I1))],[S1f(I1)./T1(I1); T1(I1)./T1(I1)],'g');
hold off

xlabel('x')
ylabel('t')

end

