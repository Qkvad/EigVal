function [y,ro]=inverzna_iteracije(A,x,sigma,tol) 
  y=x/norm(x,2);
  ro=(y'*inv(A)*y)/(y'*y);
  r=norm(inv(A)*y-ro*y,2);
  n=max(size(A));
  I=eye(n);
  while(r>tol)
    x=inv(A-sigma*I) *y;
    y=x/norm(x,2);
    ro=(y'*inv(A)*y)/(y'*y);
    r=norm(inv(A)*y-ro*y,2);
  end
end
