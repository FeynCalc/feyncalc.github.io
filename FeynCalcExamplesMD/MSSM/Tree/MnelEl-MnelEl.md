---
title: Neutralino-electron scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Mnel El -> Mnel El, MSSM, matrix element squared, tree"; 
If[$FrontEnd === Null, $FeynCalcStartupMessages = False; 
      Print[description]; ]; 
If[$Notebooks === False, $FeynCalcStartupMessages = False]; 
$LoadAddOns = {"FeynArts"}; 
Get["FeynCalc`"]
$FAVerbose = 0; 
FCCheckVersion[9, 3, 0]; 
```

![0qnnh03rto7wq](img/0qnnh03rto7wq.svg)

![02tqcun616cas](img/02tqcun616cas.svg)

![0j973yme4iv1e](img/0j973yme4iv1e.svg)

![1gj07ff4c9vo9](img/1gj07ff4c9vo9.svg)

![0yl3w9146i37j](img/0yl3w9146i37j.svg)

![173evn30flup4](img/173evn30flup4.svg)

![1qo4z5not0lhy](img/1qo4z5not0lhy.svg)

![0liutpchexhmt](img/0liutpchexhmt.svg)

![145baygm4jppw](img/145baygm4jppw.svg)

## Generate Feynman diagrams

Nicer typesetting

```mathematica
MakeBoxes[p1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(1\)]\)"; 
MakeBoxes[p2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(2\)]\)"; 
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {F[11, {1}], F[2, {1}]} -> {F[11, {1}], F[2, {1}]}, 
       InsertionLevel -> {Classes}, Model -> MSSM, 
       ExcludeParticles -> {S[1], S[2], S[3], S[4], V[_]}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1i3ojoj7w3me7](img/1i3ojoj7w3me7.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
       IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
       ChangeDimension -> 4, List -> False, SMP -> True, 
       DropSumOver -> True] //. {USf[args1__][args2__] :> 
         USf[args2, args1], Index[Sfermion, 5] :> Sfe5, 
       Conjugate[ZNeu[a__]] :> ZNeuC[a], 
       Conjugate[USf[a_, b_, c_, d_]] :> USfC[a, b, c, d]}
```

![02szolxiwy862](img/02szolxiwy862.svg)

```mathematica
Cases2[amp[0], USf]
```

![08hdphz18l07m](img/08hdphz18l07m.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, MNeu[1], SMP["m_e"], 
     MNeu[1], SMP["m_e"]]; 
```

## Evaluate the amplitude

```mathematica
amp[1] = DiracSimplify[amp[0]]; 
```

```mathematica
ampCC[1] = ComplexConjugate[amp[1], Conjugate -> 
           {ZNeuC, ZNeu, USf, USfC}] //. 
       {Conjugate[USf][a_, b_, c_, d_] :> USfC[a, b, c, d], 
         Conjugate[ZNeu][a__] :> ZNeuC[a], Conjugate[ZNeuC][a__] :> 
           ZNeu[a], Conjugate[USfC][a_, b_, c_, d_] :> 
           USf[a, b, c, d], Sfe5 -> Sfe5c}; 
```

## Square the amplitude

To avoid having too many terms, we isolate everything that is not required to calculate the spin sums

```mathematica
amp[2] = Collect2[amp[1], Spinor, LorentzIndex, 
       IsolateNames -> KK]; 
ampCC[2] = Collect2[ampCC[1], Spinor, LorentzIndex, 
       IsolateNames -> KK]; 
```

```mathematica
ampSquared[0] = DiracSimplify[FermionSpinSum[amp[2]*ampCC[2]]]; 
```

For simplicity, we neglect the masses of the external particles.

```mathematica
ampSquared[1] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
       Factor2[(#1 /. {MNeu[1] -> 0, SMP["m_e"] -> 0} & )[
           PropagatorDenominatorExplicit[FRH[ampSquared[0]]]]]]; 
```

As we will show below, the pieces that contain a Levi-Civita tensor can be discarded, so we ignore
them in the final result

```mathematica
ampSquared[2] = Simplify[SelectFree2[ampSquared[1], Eps]]
```

![1ohj9axxfy599](img/1ohj9axxfy599.svg)

The explicit dependence on the Levi-Civita tensor vanishes once we exploit the unitarity of the
sfermion mixing matrix USf

```mathematica
discarded = Simplify[SelectNotFree2[ampSquared[1], Eps]]
```

![1q9y5mah6tr7o](img/1q9y5mah6tr7o.svg)

```mathematica
Simplify[(#1 //. {USf[1, 1, re__]*USfC[1, 2, re__] :> 
              (-USf[2, 1, re])*USfC[2, 2, re], 
            USf[1, 2, re__]*USfC[1, 1, re__] :> (-USf[2, 2, re])*
                USfC[2, 1, re]} & )[Simplify[Sum[discarded, 
         {Sfe5, 1, 2}, {Sfe5c, 1, 2}]]]]
```

![1b3x5u5rst8v1](img/1b3x5u5rst8v1.svg)

## Check the final results

```mathematica
knownResults = 
     {-(SMP["e"]^4*(SMP["sin_W"]^4*(USf[Sfe5, 1, 2, 1]*
                        ((s*u*(2*s*u + t*MSf[Sfe5c, 2, 1]^2) + 
                                MSf[Sfe5, 2, 1]^2*(s*t*u + (s^2 + u^2)*
                                     MSf[Sfe5c, 2, 1]^2))*
                USf[Sfe5c, 1, 2, 1]*
                             USfC[Sfe5, 1, 2, 1] + 
               4*(2*s*u + t*MSf[Sfe5, 2, 1]^
                                    2)*(2*s*u + 
                  t*MSf[Sfe5c, 2, 1]^2)*
                             USf[Sfe5c, 2, 2, 1]*
                USfC[Sfe5, 2, 2, 1])*
                        USfC[Sfe5c, 1, 2, 1] + 4*USf[Sfe5, 2, 2, 1]*
                        ((2*s*u + t*MSf[Sfe5, 2, 1]^2)*(2*s*u + 
                                t*MSf[Sfe5c, 2, 1]^2)*
                USf[Sfe5c, 1, 2, 1]*
                             USfC[Sfe5, 1, 2, 1] + 4*(s*u*(2*s*u + 
                                    t*MSf[Sfe5c, 2, 1]^2) + 
                  MSf[Sfe5, 2, 1]^2*(
                                    
                    s*t*u + (s^2 + u^2)*MSf[Sfe5c, 2, 1]^2))*
                             USf[Sfe5c, 2, 2, 1]*
                USfC[Sfe5, 2, 2, 1])*
                        USfC[Sfe5c, 2, 2, 1])*ZNeu[1, 1]^2*
          ZNeuC[1, 1]^2 + 
                 (s*u*(2*s*u + t*MSf[Sfe5c, 2, 1]^2) + 
            MSf[Sfe5, 2, 1]^2*
                        (s*t*u + (s^2 + u^2)*MSf[Sfe5c, 2, 1]^2))*
                   SMP["cos_W"]^4*USf[Sfe5, 1, 2, 1]*
          USf[Sfe5c, 1, 2, 1]*
                   USfC[Sfe5, 1, 2, 1]*USfC[Sfe5c, 1, 2, 1]*
          ZNeu[1, 2]^2*
                   ZNeuC[1, 2]^2 + 2*SMP["cos_W"]*SMP["sin_W"]^3*
                   (USf[Sfe5, 1, 2, 
              1]*((s*u*(2*s*u + t*MSf[Sfe5c, 2, 1]^
                                      2) + 
                  MSf[Sfe5, 2, 1]^2*(s*t*u + (s^2 + u^2)*
                                     MSf[Sfe5c, 2, 1]^2))*
                USf[Sfe5c, 1, 2, 1]*
                             USfC[Sfe5, 1, 2, 1] + 
               2*(2*s*u + t*MSf[Sfe5, 2, 1]^
                                    2)*(2*s*u + 
                  t*MSf[Sfe5c, 2, 1]^2)*
                             USf[Sfe5c, 2, 2, 1]*
                USfC[Sfe5, 2, 2, 1])*
                        USfC[Sfe5c, 1, 2, 1] + 
                      2*(2*s*u + t*MSf[Sfe5, 2, 1]^2)*(2*s*u + 
                           t*MSf[Sfe5c, 2, 1]^2)*USf[Sfe5, 2, 2, 1]*
                        USf[Sfe5c, 1, 2, 1]*USfC[Sfe5, 1, 2, 1]*
                        USfC[Sfe5c, 2, 2, 1])*ZNeu[1, 1]*ZNeuC[1, 1]*
                   (ZNeu[1, 2]*ZNeuC[1, 1] + ZNeu[1, 1]*ZNeuC[1, 2]) + 
                 2*(s*u*(2*s*u + t*MSf[Sfe5c, 2, 1]^2) + 
                      MSf[Sfe5, 2, 1]^2*(s*t*u + (s^2 + u^2)*
                             MSf[Sfe5c, 2, 1]^2))*SMP["cos_W"]^3*
          SMP["sin_W"]*
                   USf[Sfe5, 1, 2, 1]*USf[Sfe5c, 1, 2, 1]*
                   USfC[Sfe5, 1, 2, 1]*USfC[Sfe5c, 1, 2, 1]*ZNeu[1, 2]*
                   ZNeuC[1, 2]*(ZNeu[1, 2]*ZNeuC[1, 1] + 
                      ZNeu[1, 1]*ZNeuC[1, 2]) + SMP["cos_W"]^2*
                   SMP["sin_W"]^2*(4*(2*s*u + t*MSf[Sfe5, 2, 1]^2)*
                        (2*s*u + t*MSf[Sfe5c, 2, 1]^2)*
             USf[Sfe5, 2, 2, 1]*
                        USf[Sfe5c, 1, 2, 1]*USfC[Sfe5, 1, 2, 1]*
                        USfC[Sfe5c, 2, 2, 1]*ZNeu[1, 1]*ZNeu[1, 2]*
                        ZNeuC[1, 1]*ZNeuC[1, 2] + USf[Sfe5, 1, 2, 1]*
                        
             USfC[Sfe5c, 1, 2, 1]*(4*(2*s*u + t*MSf[Sfe5, 2, 1]^
                                    2)*(2*s*u + 
                  t*MSf[Sfe5c, 2, 1]^2)*
                             USf[Sfe5c, 2, 2, 1]*USfC[Sfe5, 2, 2, 1]*
                ZNeu[1, 1]*
                             ZNeu[1, 2]*ZNeuC[1, 1]*ZNeuC[1, 2] + 
                           (s*u*(2*s*u + t*MSf[Sfe5c, 2, 1]^2) + 
                                MSf[Sfe5, 2, 1]^2*(s*t*u + (s^2 + u^2)*
                                     MSf[Sfe5c, 2, 1]^2))*
                USf[Sfe5c, 1, 2, 1]*
                             
                USfC[Sfe5, 1, 2, 1]*(ZNeu[1, 2]^2*ZNeuC[1, 1]^2 + 
                                
                  4*ZNeu[1, 1]*ZNeu[1, 2]*ZNeuC[1, 1]*ZNeuC[1, 2] + 
                                ZNeu[1, 1]^2*ZNeuC[1, 2]^2)))))/
         (4*(s - MSf[Sfe5, 2, 1]^2)*(-u + MSf[Sfe5, 2, 1]^2)*
            (s - MSf[Sfe5c, 2, 1]^2)*(u - MSf[Sfe5c, 2, 1]^2)*
            SMP["cos_W"]^4*SMP["sin_W"]^4)}; 
FCCompareResults[{ampSquared[2]}, knownResults, 
     Text -> {"\tCompare to FormCalc:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1mh90mfvgk0xz](img/1mh90mfvgk0xz.svg)

![1fuytza1pcbrx](img/1fuytza1pcbrx.svg)