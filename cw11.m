close all
spiking_model_660031764
%%
f=@(y)rhs(y(1:2),y(3));                                             %defining the function rhs, so it can be solved with MyTrackCurve
y0=[-90;0;-60];                                                     %setting initial conditions for MyTrackCurve
ytan=[0;0;1];
ylist=MyTrackCurve(f,0,y0,ytan,'stepsize',1);                       %running mytrackcurve to find the equilibrium points
%%
n=length(ylist);
a=zeros(1,n);                                                       %preassigning a 1xn matrix to store the types of equilibria
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);
figure(1)
hold on
for i=1:n
    if eigs(J(ylist((1:2),i),ylist(3,i)))>0
        h1=plot(ylist(3,i),ylist(1,i),'color','b','marker','o');    %classifying the type of equilibrium at each ylist point
        a(i)=1;                                                     %if it's a source, plot in blue
    elseif eigs(J(ylist((1:2),i),ylist(3,i)))<0
        h2=plot(ylist(3,i),ylist(1,i),'color','r','marker','o');
        a(i)=0;                                                     %if it's a sink, plot in red
    else
        h3=plot(ylist(3,i),ylist(1,i),'color','c','marker','o');    %if it's a saddle, plot in cyan
        a(i)=-1;
    end
end
legend([h1(1),h2(1),h3(1)],'source','sink','saddle')                %producing a legend for the plot
set(h1,'linestyle','none')
set(h2,'linestyle','none')
set(h3,'linestyle','none')

