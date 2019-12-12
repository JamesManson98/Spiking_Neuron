function df=MyJacobian(f,x,h)
    df=zeros(length(f(x)),length(x));               %preallocating a matrix df for speed
        for j=1:length(x)
            e=zeros(length(x),1);                   %cycling through matrices of length x, with a single 1 entry in the jth row.
            e(j)=1;
            df(:,j)=(1/(2*h))*(f(x+h*e)-f(x-h*e));  %using the symmetric derivative
        end
end
        
