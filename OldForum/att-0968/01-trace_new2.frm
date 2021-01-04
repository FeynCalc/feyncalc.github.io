Symbols m,u;
Vectors k,p;
AutoDeclare Index lor;
Index a,b;
Format Mathematica;

* Tr[A.GA[a].GS[k].GA[b].(1-GA[5])].Tr[B.GA[b].GS[p].GA[a].(1-GA[5])]=4*Tr[A.GS[p].(1-GA[5])]*Tr[B.GS[k].(1-GA[5])],

L A1 = g_(1,lor1);
L B1 = g_(2,lor4);

L A3 = g_(1,lor1)*g_(1,lor2)*g_(1,lor3);
L B3 = g_(2,lor4)*g_(2,lor5)*g_(2,lor6);

L LHS1 = g_(1,a)*g_(1,k)*g_(1,b)*(gi_(1)-g5_(1));
L LHS2 = g_(2,b)*g_(2,p)*g_(2,a)*(gi_(2)-g5_(2));

L RHS1 = g_(1,p)*(gi_(1)-g5_(1));
L RHS2 = g_(2,k)*(gi_(2)-g5_(2));

L res1 = A1*LHS1* B1*LHS2 - 4*A1*RHS1*B1*RHS2;

L res2 = A3*LHS1* B3*LHS2 - 4*A3*RHS1*B3*RHS2;

trace4,1;
trace4,2;
contract 0;

.sort

print res1,res2;
.end
