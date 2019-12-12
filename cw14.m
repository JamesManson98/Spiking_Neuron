close all
spiking_model_660031764
load('cw12variables')   %we load cw12variables and cw11variables
load('cw11variables')
M=@(t1,x0,p)MyIVP(@(t,x)rhs(x,p),x0,[0,t1],100);    %we define the map v
VH=ysolhopf(1); %we define the vector [VH;wH;IH] as the hopf point found in Q2
wH=ysolhopf(2);
IH=ysolhopf(3);
f=@(y)(M(y(2),[VH;y(1)],y(3))-[VH;y(1)]);   %we define the system of equations
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-6);
[~,L]=eig(J([VH;wH],IH));                   %we find the maximum imaginary eigenvalue as omega
[~,ind]=max(imag(L));
omega=imag(L(ind(1),ind(2)));
y0=[wH+1e-2;2*pi/omega;IH];                 %the time period is now found as 2*pi/omega
ytan=[1;0;0];                               %we define ytan for mytrackcurve
y=MyTrackCurve(f,0,y0,ytan);                %we find the outputs of the curve from the system of equations
J2=@(T,x,p)MyJacobian(@(x)M(T,x,p),x,1e-6); %we define the jacobian of M, used to find eigenvalues
z=zeros(1,length(y));                       %preallocating for speed
z2=zeros(1,length(y));                      %preallocating for speed
eigenvals=NaN(2,149);                       %preallocating for speed
figure(1)
hold on
for k=1:length(y) %for every trajectory, we find the eigenvalues, we can then classify the points, such that if there exists an eigenvalue greater than 1, it is unstable, otherwise it is stable
    eigenvals(:,k)=eig(J2(y(2,k),[VH;y(1,k)],y(3,k)));
    c=abs(eigenvals(:,k));
    [~,~,xt]=M(y(2,k),[VH;y(1,k)],y(3,k));  %we define the trajectories for every k
    if c(1)<=1+1e-4 && c(2)<=1+1e-4
        i1=plot(xt(1,:),xt(2,:),'color','r');   %we plot the stable trajectories red and the unstable blue
    else
        i2=plot(xt(1,:),xt(2,:),'color','b');
    end
    z(k)=max(xt(1,:));  %we also find the maximum and minimum voltage in each trajectory
    z2(k)=min(xt(1,:));
end
figure(2)
hold on
for j=1:length(y)
    eigenvals(:,j)=eig(J2(y(2,j),[VH;y(1,j)],y(3,j)));
    c=abs(eigenvals(:,j));
    aeig=[y(3,j);z(j)]; %we find the points at which the maximum and minimum voltage occurs for each current
    beig=[y(3,j);z2(j)];
    ceig=horzcat(aeig,beig);
    if c(1)<1+1e-4 && c(2)<1+1e-4
        h6=plot(ceig(1,:),ceig(2,:),'go');  %we plot stable periodic orbits green and unstable periodic orbits black
    else
        h7=plot(ceig(1,:),ceig(2,:),'ko');
    end
end
for i=1:n
    if eigs(J(ylist((1:2),i),ylist(3,i)))>0
        h1=plot(ylist(3,i),ylist(1,i),'color','b','marker','o');    %classifying the type of equilibria at each ylist point, so we can plot the bifurcation diagram
        a(i)=1;
    elseif eigs(J(ylist((1:2),i),ylist(3,i)))<0
        h2=plot(ylist(3,i),ylist(1,i),'color','r','marker','o');
        a(i)=0;
    else
        h3=plot(ylist(3,i),ylist(1,i),'color','c','marker','o');
        a(i)=-1;
    end
end
resfold=horzcat(res1fold,res2fold);
h4=plot(resfold(5,:),resfold(1,:),'kx','markersize',14);     %plotting the fold equilibria from q2
h5=plot(ysolhopf(3),ysolhopf(1),'m*','markersize',14);     %plotting the hopf bifurcation from q2
hold off
figure(3)
j1=plot(y(3,:),y(2,:)); %we plot current against time period for the periodic orbits
%%
y0=[0.2878;22.7984;22.2101]; %we use an initial guess for y0 to find our initial y0in

J=@(V,I)MyJacobian(@(V)M(V(2),[VH;V(1)],I)-[VH;V(1)],V,1e-5);   %we define the system of equations
f=@(y)[M(y(2),[VH;y(1)],y(5))-[VH;y(1)];(J(y(1:2),y(5)))*y(3:4);y(3:4)'*y(3:4)-1];
df=@(y)MyJacobian(f,y,1e-5);
res=@(y)MySolve(f,y,df,1e-6,10);


[V,~]=eig(J(y0(1:2),y0(3)));        %we define the input variables Vini, by taking the eigenvector corresponding to the minimum eigenvalues
[~,ind]=min(eig(J(y0(1:2),y0(3))));
Vini=V(:,ind);
y0in=[y0(1:2);Vini;y0(3)];
res1=res(y0in);
[~,~,xt1]=M(res1(2),[VH;res1(1)],res1(5));  %we then input our solution into the function M to obtain the trajectory of the periodic orbit at res1, the fold point of the periodic orbits

figure(2)
hold on
maxV=max(xt1(1,:)); %we find the maximum and minimum voltage obtained in the orbit
minV=min(xt1(1,:));
res1V=[res1(5);maxV];   %we then plot these onto the bifurcation diagram
res1v=[res1(5);minV];
plotminmaxV=horzcat(res1V,res1v);
h9=plot(plotminmaxV(1,:),plotminmaxV(2,:),'marker','x','markersize',14,'color','[0.5,0.5,1]');


figure(1)
hold on
i4=plot(xt1(1,:),xt1(2,:),'color','k'); %we also plot the trajectory of the orbit to compare with other trajectories in the (V,W) plane

figure(3)
hold on
j2=plot(res1(5),res1(2),'marker','*','markersize',14);  %we plot the fold of the periodic orbits on the (I0,T) plane
legend([j1,j2],'I0 vs. Period','Fold of Periodic Orbits')
res1=round(res1,4,'significant');    %rounding res1 to 4 significant figures
%%
I0=3.6; %we guess an initial I0 from the bifurcation diagram, close to the saddle connection
f5=@(ysaddle)rhs(ysaddle(1:2),ysaddle(3));  %we define the system of equations to find the equilibrium4_2 at a fixed I0 (I0=3.6) the equilibrium to be used in our initial guess
f6=@(ysaddle)(I0-ysaddle(3));
feq2=@(y)[f5(y);f6(y)];
dfeq2=@(y)MyJacobian(feq2,y,1e-5);
for i=1:length(ylist)-1
    if ylist(3,i+1)<=I0 && ylist(3,i)>I0 && a(i)==-1
        yguess2=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0];   %we guess a value using ylist found in Q1
    end
end
equilibrium4_2=MySolve(feq2,yguess2,dfeq2,1e-6,10); %we have now found the separatrix for I0=3.6 which is close to p, the value of I for the connecting orbit.
M=@(t1,x0,p)MyIVP(@(t,x)rhs(x,p),x0,[0,t1],600);
Ts=0.93*(35);   %we define Ts as the modulues of the time taken to reach the greatest x1 on the trajectory going backwards from the stable separatrix
Tu=0.99*(100);  %we define Tu as the time taken to reach the greatest x1 on the trajectory going forwards to the unstable separatrix
x2=0.22;        %we find x2, an estimate of the value at which the greatest x1 on the trajectory occurs
s=0.1;          %we find a suitable s, to move from the separatrix
select=@(M,r,c)M(r,c);  %we use the select function to create a function that outputs only the top or bottom row of M
M2=@(t1,x0,p)select(M(t1,x0,p),2,1);    %we define M2 as the bottom row of the output of M
M1=@(t1,x0,p)select(M(t1,x0,p),1,1);    %we define M1 as the top row of the output of M
f=@(y)rhs(y(1:2),y(5));                 %we define the system of equations feq, in terms of y              
f1=@(y)M2(-y(3),[y(1);y(2)]+s*vsvalue([y(1);y(2);y(5)],rhs),y(5))-x2;
f2=@(y)M2(y(4),[y(1);y(2)]+s*vuvalue([y(1);y(2);y(5)],rhs),y(5))-x2;
f3=@(y)M1(-y(3),[y(1);y(2)]+s*vsvalue([y(1);y(2);y(5)],rhs),y(5))-M1(y(4),[y(1);y(2)]+s*vuvalue([y(1);y(2);y(5)],rhs),y(5));
feq=@(y)[f(y);f1(y);f2(y);f3(y)];
dfeq=@(y)MyJacobian(feq,y,1e-6);        %we define dfeq as the Jacobian of the system of equations, to be input in MySolve
y0_saddle=[equilibrium4_2(1);equilibrium4_2(2);Ts;Tu;I0];   %we find y0_saddle, the initial guess for y, from our values found above
yout=MySolve(feq,y0_saddle,dfeq,1e-6,10);   %we run MySolve and obtain an output the matrix yout, the solution
[~,~,xt_saddle]=M(200,[yout(1);yout(2)]+s*vuvalue([yout(1);yout(2);yout(5)],rhs),yout(5));  %we then find the trajectory of the periodic orbit at yout
z3=max(xt_saddle(1,:)); %we find the maximum and minimum voltage
z4=min(xt_saddle(1,:));
figure(2)
hold on
h8=plot(yout(5),z3,'marker','x','color','r');   %we plot the maximum and minimum voltage on the bifurcation diagram
plot(yout(5),z4,'marker','x','color','r')
legend([h1,h2,h3,h4,h5,h6,h7,h9,h8],'source','sink','saddle','fold equilibrium','hopf','stable','unstable','fold in periodic orbit','connecting orbit') %making the legend
figure(1)
hold on
i3=plot(xt_saddle(1,:),xt_saddle(2,:),'color','g'); %plotting the connecting orbits on the graph of phase portraits for all periodic orbits
legend([i1,i2,i4,i3],'stable','unstable','fold of periodic orbits','connecting orbit') %making the legend
yout=round(yout,4,'significant');
