Symbols m,u;
Vectors k1,k2,p,p1,p2,q,s;
AutoDeclare Index lor;
Format Mathematica;

L TrA2B1 = ((-(m*gi_(1))+g_(1,p2))*g_(1,lor4)*(m*gi_(1)+g_(1,p1))*g_(1,lor3)*(-(m*gi_(1))+g_(1,p))*g_(1,lor2)*g_(1,k2)*g_(1,lor1)*
(-g5_(1)+gi_(1))*g_(2,k1)*g_(2,lor2)*(-(u*gi_(2))-g_(2,p1)-g_(2,p2)+g_(2,q))*g_(2,lor4)*(-(u*gi_(2))+g_(2,q))*
(gi_(2)+g5_(2)*g_(2,s))*g_(2,lor3)*(-(u*gi_(2))-g_(2,p)-g_(2,p1)+g_(2,q))*g_(2,lor1)*(-g5_(2)+gi_(2)));

L TrA2B2 = ((-(m*gi_(1))+g_(1,p))*g_(1,lor4)*(m*gi_(1)+g_(1,p1))*g_(1,lor3)*(-(m*gi_(1))+g_(1,p2))*g_(1,lor2)*g_(1,k2)*g_(1,lor1)*
(-g5_(1)+gi_(1))*g_(2,k1)*g_(2,lor2)*(-(u*gi_(2))-g_(2,p)-g_(2,p1)+g_(2,q))*g_(2,lor4)*(-(u*gi_(2))+g_(2,q))*
(gi_(2)+g5_(2)*g_(2,s))*g_(2,lor3)*(-(u*gi_(2))-g_(2,p1)-g_(2,p2)+g_(2,q))*g_(2,lor1)*(-g5_(2)+gi_(2)));

L resFL = TrA2B1 + TrA2B2;

trace4,1;
trace4,2;
contract 0;

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
