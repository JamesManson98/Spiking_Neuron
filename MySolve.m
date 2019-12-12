function [x,converged,J]=MySolve(f,x0,df,tol,maxit)
    xold=x0;                                      
    converged=0;                                    %originally the function has not converged
    for i=1:maxit
        x=xold-df(xold)\f(xold);                    %we compute a new x via Newton iteration
        if norm(x-xold)<tol && norm(f(xold))<tol    %if the correction and the residual are less than tol, it has converged and we break
            converged=1;
            break
        end
        xold=x;                                     %we set x as xold to find the next value of x
    end
    J=df(xold);                                     %we find the last Jacobian after the for loop and output it
end