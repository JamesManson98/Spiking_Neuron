%%
close all  %closing all previously plotted graphs
spiking_model_660031764
load('cw11variables')   %loading variables from cw11 so ylist entries can be used
%finding the initial guesses (y1fold,y2fold,yhopf), where the type of equilibria changes
for i=1:n-1
    if a(i)==0 && a(i+1)==-1
        y1fold=ylist(:,i);
    end
    if a(i)==-1 && a(i+1)==1
        y2fold=ylist(:,i);
    end
    if a(i)==1 && a(i+1)==0
        yhopf=ylist(:,i);
    end
end
%defining res in terms of functions f and df using MyJacobian and MySolve
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);
f=@(y)[rhs(y(1:2),y(5));(J(y(1:2),y(5)))*y(3:4);y(3:4)'*y(3:4)-1];
df=@(y)MyJacobian(f,y,1e-5);
res=@(y)MySolve(f,y,df,1e-6,10);

%finding the first fold
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);                %defining the input variable y1foldin
[V,~]=eig(J(y1fold(1:2),y1fold(3)));
[~,ind]=min(eig(J(y1fold(1:2),y1fold(3))));
Vini=V(:,ind);                                          %finding the initial vector Vini, corresponding to the minimum eigenvalue
y1foldin=[y1fold(1:2);Vini;y1fold(3)];
res1fold=res(y1foldin);                                 %finding the output for the first fold
%finding the second fold
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);                %defining the input variable y2foldin
[V,L]=eig(J(y2fold(1:2),y2fold(3)));
[~,ind]=min(eig(J(y2fold(1:2),y2fold(3))));
Vini=V(:,ind);                                          %finding the initial vector Vini, corresponding to the minimum eigenvalue
y2foldin=[y2fold(1:2);Vini;y2fold(3)];
res2fold=res(y2foldin);                                 %finding the output for the second fold
%%
%finding the hopf bifurcation
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);                %setting up the system of equations with f and df to be solved by MySolve
f=@(y)[rhs(y(1:2),y(3));trace(J(y(1:2),y(3)))];
df=@(y)MyJacobian(f,y,1e-5);
ysolhopf=MySolve(f,yhopf,df,1e-6,10);
hold on
for i=1:n
    if a(i)==1
        h1=plot(ylist(3,i),ylist(1,i),'color','b','marker','o');    %classifying the type of equilibrium at each ylist point
    elseif a(i)==0
        h2=plot(ylist(3,i),ylist(1,i),'color','r','marker','o');
    else
        h3=plot(ylist(3,i),ylist(1,i),'color','c','marker','o');
    end
end
h4=plot(res1fold(5),res1fold(1),'kx','markersize',14);     %plotting the first fold
h6=plot(ysolhopf(3),ysolhopf(1),'m*','markersize',14);     %plotting the hopf bifurcation
h5=plot(res2fold(5),res2fold(1),'kx','markersize',14);     %plotting the second fold
legend([h1,h2,h3,h4,h6],'source','sink','saddle','fold','hopf') %making the legend
hold off
res1foldround=round(res1fold,4,'significant');
res2foldround=round(res2fold,4,'significant');
ysolhopfround=round(ysolhopf,4,'significant');

