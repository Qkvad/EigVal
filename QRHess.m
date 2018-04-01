function [A,pomA]=QRHess(A)
  A=hessenberg(A);
  pomA=A;
  n=max(size(A));
  Q=eye(n);
  c=A(1,1)/sqrt(A(1,1)^2 + A(2,1)^2);
  s=A(2,1)/sqrt(A(1,1)^2 + A(2,1)^2);
  Q(1,1)=c; Q(2,2)=c; Q(1,2)=s; Q(2,1)=-s;
  A=Q*A*Q';
  for j=1:1000
    for i=1:n-2
      if A(i+1,i)~=0
        c=A(i+1,i)/sqrt(A(i+1,i)^2 + A(i+2,i)^2);
        s=A(i+2,i)/sqrt(A(i+1,i)^2 + A(i+2,i)^2);
      else
        c=0; s=1;
      endif
      Q=eye(n);
      Q(i+1,i+1)=c; Q(i+2,i+2)=c; Q(i+1,i+2)=s; Q(i+2,i+1)=-s;
      A=Q*A*Q';
    end
  end
end
