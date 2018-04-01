function [y,ro]=test_II(A,x,sigma,tol)
  n=max(size(A));
  t=sort(rand(n,1),'descend');
  D=diag(t);
  [Q,R]=qr(A);
  A=Q'*D*Q;
  [y,ro]=inverzna_iteracije(A,x,sigma,tol);
end
