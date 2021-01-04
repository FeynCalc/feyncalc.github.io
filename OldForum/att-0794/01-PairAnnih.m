(* ::Package:: *)

$LoadPhi=False;
$LoadFeynArts=True;
<<HighEnergyPhysics`FeynCalc`

$FAVerbose=0;
tops=CreateTopologies[0,2->2];
Paint[tops,AutoEdit->False,ColumnsXRows->{4,1}];

inserttops=InsertFields[tops,{F[2,{1}],-F[2,{1}]}->{V[1],V[1]},InsertionLevel->{Classes}];
Paint[inserttops,AutoEdit->False,ColumnsXRows->{3,1}];

Amp=Map[ReplaceAll[#,FeynAmp[_,_,amp_,___]:>amp]&,Apply[List,CreateFeynAmp[inserttops,Truncated->False]]]//.{(a1__ DiracGamma[6] a2__+a1__ DiracGamma[7] a2__):>a1 a2,NonCommutative[x___]->x,FermionChain->DOT,FourMomentum[Outgoing,1]->k1,FourMomentum[Outgoing,2]->k2,FourMomentum[Incoming,1]->p1,FourMomentum[Incoming,2]->p2,DiracSpinor->Spinor,LorentzIndex[Index[Lorentz,x_]]:>LorentzIndex[ToExpression["Lor"<>ToString[x]]],Conjugate[PolarizationVector][_,x_,y_]:>Conjugate[PolarizationVector[x,y]]}

SetMandelstam[s,t,u,p1,p2,-k1,-k2,ME,ME,0,0];

SqAmp=Total[Amp] Total[(ComplexConjugate[Amp]//FCRenameDummyIndices)]//PropagatorDenominatorExplicit//Expand//ReplaceRepeated[#,{Pair[LorentzIndex[x1_],Momentum[Polarization[y_,I]]]Pair[LorentzIndex[x2_],Momentum[Polarization[y_,-I]]]
   :>-MT[x1,x2]}]&//Contract//FermionSpinSum[#,ExtraFactor->1/2^2,SpinorCollect->True]&//ReplaceAll[#,DiracTrace[x___]:>DiracTrace[x,DiracTraceEvaluate->True]]&

ClearScalarProducts
ScalarProduct[k1,k1]=0;
ScalarProduct[k2,k2]=0;
ScalarProduct[p1,p1]=ME^2;
ScalarProduct[p2,p2]=ME^2;

SqAmp2=((SqAmp/.{s->-t-u+2ME^2})/.{t->ScalarProduct[p1-k1],u->ScalarProduct[p1-k2]})//ExpandScalarProduct//Simplify//Expand


2EL^4 Collect[Expand[SqAmp2/(2EL^4)],{ME^2,ME^4}]

