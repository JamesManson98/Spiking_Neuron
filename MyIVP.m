function [xend,t,xt]=MyIVP(f,x0,tspan,N)
    b=size(x0);                                 %we find the size of x0
    c=b(1);                                     %we find the length of x0
    h=(tspan(2)-tspan(1))/N;                    %we find the increment to increase t by each time
    t=tspan(1):h:tspan(2);                      %creating a matrix with t ranging between tspan(1) and tspan(2) in increment h
    xt=zeros(c,N+1);                            %preallocating for speed
    xt(:,1)=x0;
    for i=1:N
        k1=h*f(t(i),xt(:,i));                   %we follow the Dormand-Prince method to explicitly solve the ODE for every i=1:N
        k2=h*f(t(i)+(1/5)*h,xt(:,i)+(1/5)*k1);
        k3=h*f(t(i)+(3/10)*h,xt(:,i)+(3/40)*k1+(9/40)*k2);
        k4=h*f(t(i)+(4/5)*h,xt(:,i)+(44/45)*k1-(56/15)*k2+(32/9)*k3);
        k5=h*f(t(i)+(8/9)*h,xt(:,i)+(19372/6561)*k1-(25360/2187)*k2+(64448/6561)*k3-(212/729)*k4);
        k6=h*f(t(i)+h,xt(:,i)+(9017/3168)*k1-(355/33)*k2+(46732/5247)*k3+(49/176)*k4-(5103/18636)*k5);
        xt(:,i+1)=xt(:,i)+(35/384)*k1+(500/1113)*k3+(125/192)*k4-(2187/6784)*k5+(11/84)*k6;     %we track the curve by inputting every row into xt
    end
    xend=xt(:,i+1);                             %we find the last value of xt and output as xend
end