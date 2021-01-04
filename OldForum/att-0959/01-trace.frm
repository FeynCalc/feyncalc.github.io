Symbol u,m;
Vectors k1,p,p1,p2,q,s,k2;
Index alpha,alpha1,beta,beta1;
Format Mathematica;

L Line13 = ((-(m*gi_(1))+g_(1,p2))*g_(1,beta1)*(m*gi_(1)+g_(1,p1))*g_(1,beta)*(-(m*gi_(1))+g_(1,p))*g_(1,alpha1)*g_(1,k2)*g_(1,alpha)*
(-g5_(1)+gi_(1)));

L Line14 = (g_(1,k1)*g_(1,alpha1)*(-(u*gi_(1))-g_(1,p1)-g_(1,p2)+g_(1,q))*g_(1,beta1)*(-(u*gi_(1))+g_(1,q))*(gi_(1)+g5_(1)*g_(1,s))*
g_(1,beta)*(-(u*gi_(1))-g_(1,p)-g_(1,p1)+g_(1,q))*g_(1,alpha)*(-g5_(1)+gi_(1)));

L Line15 = ((-(m*gi_(1))+g_(1,p))*g_(1,beta1)*(m*gi_(1)+g_(1,p1))*g_(1,beta)*(-(m*gi_(1))+g_(1,p2))*g_(1,alpha1)*g_(1,k2)*g_(1,alpha)*
(-g5_(1)+gi_(1)));

L Line16 = (g_(1,k1)*g_(1,alpha1)*(-(u*gi_(1))-g_(1,p)-g_(1,p1)+g_(1,q))*g_(1,beta1)*(-(u*gi_(1))+g_(1,q))*(gi_(1)+g5_(1)*g_(1,s))*
g_(1,beta)*(-(u*gi_(1))-g_(1,p1)-g_(1,p2)+g_(1,q))*g_(1,alpha)*(-g5_(1)+gi_(1)));

L resFL = Line13*Line14+Line15*Line16;
trace4,1;
contract;
repeat;
  id p.p = m^2;
  id p1.p1 = m^2;
  id p2.p2 = m^2;
  id k1.k1 = 0;
  id k2.k2 = 0;
  id q.q = u^2;
  id q.s = 0;
endrepeat;

.sort
print resFL;
.end
