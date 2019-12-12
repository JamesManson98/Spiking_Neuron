function vs=vsvalue(y,rhs)
J=@(x,I0)MyJacobian(@(x)rhs(x,I0),x,1e-5);  
[V,~]=eig(J(y(1:2),y(3)));                  %finding eigenvectors corresponding to the maximum eigenvalues
[~,ind]=min(real(eig(J(y(1:2),y(3)))));
vs=V(:,ind);    
end