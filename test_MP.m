function [y,ro]=test_MP(A,x,tol)
  n=max(size(A));
  t=sort(rand(n,1),'descend');
  D=diag(t);
  [Q,R]=qr(A);
  A=Q*D*Q';
  [y,ro]=metoda_potencija(A,x,tol);
end
