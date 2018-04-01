function [X]=iteracije_potprostora(A,X,l)
  [Q,R]=qr(X);
  k=0;
  while(k<1000)
    X=A*Q(:,1:l);
    [Q,R]=qr(X);
    k=k+1;
  end
end
