function [A,U]=hessenberg(A)
  n=length(A(:,1));
  U=eye(n);
  for i=1:n-2
    v=A(i+1:n,i);
    if(v(1)>0)
      v(1)=A(i+1,i)+norm(A(i+1:n,i),2);
    else
      v(1)=v(1)-norm(A(i+1:n,i),2);
    end
    H=eye(n);
    H(i+1:n,i+1:n)=eye(n-i)-2*v*v'/(v'*v);
    U=H*U;
    A=H*A*H';
  end
  for i=2:n
    for j=1:i
      if(abs(A(i,j))<0.000001) 
        A(i,j)=0;
      end
    end
  end   
end
