close all
spiking_model_660031764
load('cw11variables')
%% finding initial guesses
%the following finds initial guesses by cycling through the length of
%ylist and finding the points where the stability changes, and taking the
%mean of the two points either side.
for i=1:n-1
    I0_1=-50;
    I0_2=-30;
    I0_3=15;
    I0_4=30;
    I0_5=5;
    I0_6=22;
    if ylist(3,i)<=I0_1 && ylist(3,i+1)>I0_1
        yguess1=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_1];
    elseif ylist(3,i)<=I0_2 && ylist(3,i+1)>I0_2 && a(i)==0
        yguess2=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_2];
    elseif ylist(3,i+1)<=I0_2 && ylist(3,i)>I0_2 && a(i)==-1
        yguess3=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_2];
    elseif ylist(3,i)<=I0_2 && ylist(3,i+1)>I0_2 && a(i)==1
        yguess4=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_2];
    elseif ylist(3,i)<=I0_3 && ylist(3,i+1)>I0_3
        yguess5=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_3];
    elseif ylist(3,i)<=I0_4 && ylist(3,i+1)>I0_4
        yguess6=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_4];
    elseif ylist(3,i)<=I0_5 && ylist(3,i+1)>I0_5 && a(i)==0
        yguess7=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_5];
    elseif ylist(3,i+1)<=I0_5 && ylist(3,i)>I0_5 && a(i)==-1
        yguess8=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_5];
    elseif ylist(3,i)<=I0_5 && ylist(3,i+1)>I0_5 && a(i)==1
        yguess9=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_5];
    elseif ylist(3,i)<=I0_6 && ylist(3,i+1)>I0_6
        yguess10=[(ylist(1,i)+ylist(1,i+1))/2;(ylist(2,i)+ylist(2,i+1))/2;I0_6];
    end
end
%% solving to find equilibria
feq=@(y)rhs(y(1:2),y(3));
dfeq=@(y)MyJacobian(feq,y,1e-5);
equilibrium1=MySolve(feq,yguess1,dfeq,1e-6,10); %sink
equilibrium2=MySolve(feq,yguess2,dfeq,1e-6,10); %sink
equilibrium3=MySolve(feq,yguess3,dfeq,1e-6,10); %saddle
equilibrium4=MySolve(feq,yguess4,dfeq,1e-6,10); %source
equilibrium5=MySolve(feq,yguess5,dfeq,1e-6,10); %source
equilibrium6=MySolve(feq,yguess6,dfeq,1e-6,10); %sink
f1=@(y)rhs(y(1:2),y(3));    %we implement additional functions to fix I0, as otherwise our guess for a point at I0 will be inaccurate
f2=@(y)(I0_5-y(3));
feq=@(y)[f1(y);f2(y)];
dfeq=@(y)MyJacobian(feq,y,1e-5);
equilibrium7=MySolve(feq,yguess7,dfeq,1e-6,10); %sink
equilibrium8=MySolve(feq,yguess8,dfeq,1e-6,10); %saddle
equilibrium9=MySolve(feq,yguess9,dfeq,1e-6,10); %source
f2=@(y)(I0_6-y(3));
feq=@(y)[f1(y);f2(y)];
equilibrium10=MySolve(feq,yguess10,dfeq,1e-6,10);   %sink
%% example of a phase portrait that occurs at I0=-50
%sink at I0=-50
figure(1)
hold on
I0=I0_1;
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK))); %implementing the equation for the V nullcline
V=-91:0.1:40;
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));   %defining the V nullcline and W nullcline
end
h1=plot(V(:),W(:));     %plotting the V nullcline
h2=plot(V(:),W2(:));    %plotting the W nullcline
h3=plot(equilibrium1(1),equilibrium1(2),'marker','*','markersize',14);  %plotting the sink
V0=[yguess1(1);yguess1(2)];         %implementing the guess
[~,~,xt1]=MyIVP(@(t,x)rhs(x,I0),V0,[0,-100],300);   %solving to find a trajectory in MyIVP
h4=plot(xt1(1,:),xt1(2,:),'.-');    %plotting the trajectory
hold off
xlim([-91,40])
ylim([-0.2,1])
legend([h1,h2,h3,h4],'V nullcline','W nullcline','Sink','Trajectory')
%% example of a phase portrait that occurs at I0=-30
%Nullclines for I0=-30
I0=I0_2;
figure(2)
hold on
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK))); %implementing the nullclines as above
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));
end
h1=plot(V(:),W(:));
h2=plot(V(:),W2(:));
xlim([-91,40])
ylim([-0.2,1])

%using the initial guess and MyIVP to find a sink for I0=-30
V0=[yguess2(1);yguess2(2)];
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,-100],300);
h3=plot(xt(1,:),xt(2,:),'.-');

%using the initial guess and MyIVP to find a source for I0=-30
V0=[yguess4(1);yguess4(2)];
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,1000],30000);
h4=plot(xt(1,:),xt(2,:),'.-');

%saddle for I0=-30
y=equilibrium3;         %we use the equilibrium to define a guess for where to begin the trajectories
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);
[V,~]=eig(J(y(1:2),y(3)));
[min1,ind]=min(eig(J(y(1:2),y(3))));
[max1,ind2]=max(eig(J(y(1:2),y(3))));
V1=V(:,ind);            %we define the eigenvectors corresponding to the minimum and maximum eigenvalues
V2=V(:,ind2);
y0=y(1:2)+(1e-1)*V1;    %we derive guesses for the beginning of the trajectories
y1=y(1:2)+(1e-1)*V2;
y2=y(1:2)-(1e-1)*V1;
y3=y(1:2)-(1e-1)*V2;
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y0,[0,-100],300);    %we integrate backwards for stable separatrices
h5=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y2,[0,-100],300);
h6=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y1,[0,100],300);     %we integrate forwards for unstable separatrices
h7=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y3,[0,100],300);
h8=plot(xt(1,:),xt(2,:),'.-');                      %we then plot all the trajectories and equilibria
h9=plot(equilibrium2(1),equilibrium2(2),'marker','*','markersize',14);
h10=plot(equilibrium3(1),equilibrium3(2),'marker','o','markersize',14);
h11=plot(equilibrium4(1),equilibrium4(2),'marker','+','markersize',14);
legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11],'V nullcline','W nullcline','Trajectory 1','Trajectory 2','Trajectory 3','Trajectory 4','Trajectory 5','Trajectory 6','Sink','Saddle','Source')
hold off
%% example of a phase portrait that occurs for 6.5424<I0<21.8776
%source for I0=15
figure(3)
I0=I0_3;
V0=[yguess5(1)-10;yguess5(2)];  %we find an initial guess (-10 as it converges to the PO quicker)
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,4000],120000);
h1=plot(xt(1,:),xt(2,:),'.-','color','b');
hold on
h2=plot(xt(1,119000:120001),xt(2,119000:120001),'.-','color','g'); %we plot the last values as an estimate of the periodic orbit
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK))); %we implement the nullclines as above
V=-91:0.1:40;
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));
end
h3=plot(V(:),W(:));
h4=plot(V(:),W2(:));
h5=plot(equilibrium5(1),equilibrium5(2),'marker','+','markersize',14);  %we plot the equilibria, nullclines and the trajectories
xlim([-91,40])
ylim([-0.2,1])
legend([h1,h2,h3,h4,h5],'Trajectory 1','Periodic Orbit','V nullcline','W nullcline','Source')
hold off
%% example of a phase portrait that occurs for I0=30
%sink for I0=30
figure(4)
I0=I0_4;
V0=[yguess6(1);yguess6(2)]; %we use an initial guess and integrate backwards away from the sink
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,-1000],30000);
h1=plot(xt(1,:),xt(2,:),'.-','color','b');
hold on
V=-91:0.1:40;   %we find the nullclines as above
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK)));
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));
end
h2=plot(V(:),W(:));
h3=plot(V(:),W2(:));
h4=plot(equilibrium6(1),equilibrium6(2),'marker','*','markersize',14); %we plot equilibria, nullclines, and trajectories
legend([h1,h2,h3,h4],'Trajectory 1','V nullcline','W nullcline','Sink')
xlim([-91,40])
ylim([-0.2,1])
hold off

%% example of a phase portrait (I0=5)
%Nullclines for I0=5
I0=I0_5;
figure(5)
hold on
V=-91:0.001:40; %we find the nullclines as above
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK)));
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));
end
h1=plot(V(:),W(:));
h2=plot(V(:),W2(:));
xlim([-91,40])
ylim([-0.2,1])

%we find a sink for I0=5 by integrating backwards with MyIVP
V0=[yguess7(1);yguess7(2)];
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,-100],300);
h3=plot(xt(1,:),xt(2,:),'.-');

%we find a source for I0=5 by integrating forwards with MyIVP
V0=[yguess9(1);yguess9(2)];
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,1000],30000);
h4=plot(xt(1,:),xt(2,:),'.-');

%saddle for I0=5
y=equilibrium8;     %we use the equilibrium to guess where to start the trajectory
J=@(x,I)MyJacobian(@(x)rhs(x,I),x,1e-5);
[V,L]=eig(J(y(1:2),y(3)));
[min2,ind]=min(eig(J(y(1:2),y(3))));    %we define the eigenvectors corresponding to the minimum and maximum eigenvalues
[max2,ind2]=max(eig(J(y(1:2),y(3))));
V1=V(:,ind);            %we define the eigenvectors corresponding to the minimum and maximum eigenvalues
V2=V(:,ind2);
y0=y(1:2)+(1e-1)*V1;
y1=y(1:2)+(1e-1)*V2;
y2=y(1:2)-(1e-1)*V1;
y3=y(1:2)-(1e-1)*V2;
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y0,[0,-100],300);    %we integrate backwards for stable separatrices
h5=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y2,[0,-100],300);
h6=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y1,[0,100],300);     %we integrate forwards for unstable separatrices
h7=plot(xt(1,:),xt(2,:),'.-');
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),y3,[0,100],300);
h8=plot(xt(1,:),xt(2,:),'.-');
h9=plot(equilibrium7(1),equilibrium7(2),'marker','*','markersize',14);  %we then plot all trajectories and equilibria
h10=plot(equilibrium8(1),equilibrium8(2),'marker','o','markersize',14);
h11=plot(equilibrium9(1),equilibrium9(2),'marker','+','markersize',14);
legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11],'V nullcline','W nullcline','Trajectory 1','Trajectory 2','Trajectory 3','Trajectory 4','Trajectory 5','Trajectory 6','Sink','Saddle','Source')
hold off

%% example of a phase portrait I0=22 (sink)
figure(6)
I0=I0_6;
V0=[yguess10(1);yguess10(2)];       %we use the guess for y derived from ylist
[~,~,xt]=MyIVP(@(t,x)rhs(x,I0),V0,[0,-1000],30000); %we then integrate backwards (away from the sink)
h1=plot(xt(1,:),xt(2,:),'.-','color','b');
hold on
h2=plot(xt(1,29000:30001),xt(2,29000:30001),'.-','color','g');  % we plot the last values (periodic orbit) green
V=-91:0.001:40;     %we define the nullclines as above
Vnull=@(V)((I0-gL*(V-VL)-gCa*minf(V)*(V-VCa))/(gK*(V-VK)));
for i=1:length(V)
    W(i)=Vnull(V(i));
    W2(i)=winf(V(i));
end
h3=plot(V(:),W(:));
h4=plot(V(:),W2(:));
xlim([-91,40])
ylim([-0.2,1])
h5=plot(equilibrium10(1),equilibrium10(2),'marker','*','markersize',14);    %we plot the sink
xlim([-91,40])
ylim([-0.2,1])
legend([h1,h2,h3,h4,h5],'Trajectory 1','Periodic Orbit','V nullcline','W nullcline','Sink')
hold off