function ylist=MyTrackCurve(userf,~,y0,ytan,varargin)
    default={'nmax',150,'stepsize',1e-2,'smin',1e-4,'smax',2,'h',1e-4,'tol',1e-6,'maxit',10,'stop',@(y)y(end)>50};  %setting default matrix for optional inputs
    options=MySetOptions(default,varargin); %calling optional inputs
    s=0;                                    %setting the initial stepsize to 0
    smin=options.smin;                      %implementing the options
    smax=options.smax;
    h=options.h;                            
    tol=options.tol;
    maxit=options.maxit;
    yj=y0;                                  %setting the initial y as the y used within the for loop
    n=length(userf(y0));
    b=zeros(n+1,1);                         
    b(n+1)=1;                               %setting the matrix as 1 in the bottom element, with the rest 0s
    k=1;
    while k<options.nmax                    %where the steps taken is less than nmax, we take another step
        ypred=yj+s*ytan;                    %setting the predictor for y
        eq1=@(y)(ytan.'*(y-ypred));         %setting up the system of equations
        f=@(y)[userf(y);eq1(y)];
        df=@(y)MyJacobian(f,y,h);           %finding the Jacobian to be input into MySolve
        [y,converged,J]=MySolve(f,yj,df,tol,maxit); %Running mysolve to find y to output in ylist
        if converged==0                     %if it doesn't converge we vary the stepsize until it converges
            s=s/2;
            if s<smin                       %if s<smin, the minimum stepsize then we return
                return
            end
            continue
        end
        if s==0                             %if s is on the initial step, where we use stepsize 0, we set the step to stepsize/1.2 (as it's later multiplied by 1.2)
            s=options.stepsize/1.2;
        end
        s=min([s*1.2,smax]);                %we increase the stepsize after it successfully converges
        if ~options.stop(y)
            ylist(:,k)=y;                   %if we are not at the condition where y stops, then we output the calculated y in ylist
        else
            return                          %if we are at the condition where y stops, we return
        end
        z=J\b;                              %we calculate the new ytan
        znew=z/norm(z,inf);                 
        z2=sign(z.'*ytan);
        ytan=znew*z2;
        yj=y;                               %we set the new y as the output from MySolve
        k=k+1;                              %we add 1 to the number of iterations occurred
    end
end
    
    
    

