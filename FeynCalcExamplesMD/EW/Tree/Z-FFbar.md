---
title: Z boson decaying into a fermion-antifermion pair
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Z -> F Fbar, EW, total decay rate, tree"; 
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

```mathematica
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
MakeBoxes[l1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(l\), \(1\)]\)"; 
MakeBoxes[l2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(l\), \(2\)]\)"; 
MakeBoxes[p1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(1\)]\)"; 
MakeBoxes[p2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(2\)]\)"; 
```

```mathematica
diagsLeptonsMassless = InsertFields[CreateTopologies[0, 
         1 -> 2], {V[2]} -> {F[1, {l}], -F[1, {l}]}, 
       InsertionLevel -> {Classes}]; 
Paint[diagsLeptonsMassless, ColumnsXRows -> {2, 1}, 
     Numbering -> Simple, SheetHeader -> None, 
     ImageSize -> {512, 256}]; 
```

![1oygufewbrsfs](img/1oygufewbrsfs.svg)

```mathematica
diagsLeptonsMassive = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[2]} -> {F[2, {l}], -F[2, {l}]}, InsertionLevel -> 
         {Classes}]; 
Paint[diagsLeptonsMassive, ColumnsXRows -> {2, 1}, 
     Numbering -> Simple, SheetHeader -> None, 
     ImageSize -> {512, 256}]; 
```

![1q7vtdy7ifavx](img/1q7vtdy7ifavx.svg)

```mathematica
diagsUpQuarks = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[2]} -> {F[3, {l}], -F[3, {l}]}, InsertionLevel -> 
         {Classes}]; 
Paint[diagsUpQuarks, ColumnsXRows -> {2, 1}, 
     Numbering -> Simple, SheetHeader -> None, 
     ImageSize -> {512, 256}]; 
```

![00xkjkgx493x9](img/00xkjkgx493x9.svg)

```mathematica
diagsDownQuarks = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[2]} -> {F[4, {l}], -F[4, {l}]}, InsertionLevel -> 
         {Classes}]; 
Paint[diagsDownQuarks, ColumnsXRows -> {2, 1}, 
     Numbering -> Simple, SheetHeader -> None, 
     ImageSize -> {512, 256}]; 
```

![18kf6w1q2q1dn](img/18kf6w1q2q1dn.svg)

## Obtain the amplitudes

```mathematica
ampLeptonsMassless[0] = FCFAConvert[CreateFeynAmp[
       diagsLeptonsMassless], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {l1, l2}, List -> False, 
     ChangeDimension -> 4, DropSumOver -> True, SMP -> True, 
     Contract -> True, FinalSubstitutions -> 
       {MLE[l] -> SMP["m_l"], SMP["e"] -> 
           Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_Z"]^2*SMP["cos_W"]^2*
               SMP["sin_W"]^2]}]
```

![1om9p6yrdv9dj](img/1om9p6yrdv9dj.svg)

```mathematica
ampLeptonsMassive[0] = FCFAConvert[CreateFeynAmp[
       diagsLeptonsMassive], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p1, p2}, List -> False, 
     ChangeDimension -> 4, DropSumOver -> True, SMP -> True, 
     Contract -> True, FinalSubstitutions -> 
       {MLE[l] -> SMP["m_l"], SMP["e"] -> 
           Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_Z"]^2*SMP["cos_W"]^2*
               SMP["sin_W"]^2]}]
```

![0mc2busu5ypyx](img/0mc2busu5ypyx.svg)

```mathematica
ampUpQuarks[0] = FCFAConvert[CreateFeynAmp[diagsUpQuarks], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, FinalSubstitutions -> 
       {MQU[l] -> SMP["m_q"], SMP["e"] -> 
           Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_Z"]^2*SMP["cos_W"]^2*
               SMP["sin_W"]^2]}]
```

![1l0up61lhr6oc](img/1l0up61lhr6oc.svg)

```mathematica
ampDownQuarks[0] = FCFAConvert[CreateFeynAmp[diagsDownQuarks], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, FinalSubstitutions -> 
       {MQD[l] -> SMP["m_q"], SMP["e"] -> 
           Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_Z"]^2*SMP["cos_W"]^2*
               SMP["sin_W"]^2]}]
```

![1scsgnw42cnd8](img/1scsgnw42cnd8.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]
SP[p] = SMP["m_Z"]^2; 
SP[k1] = SMP["m_q"]^2; 
SP[k2] = SMP["m_q"]^2; 
SP[p1] = SMP["m_l"]^2; 
SP[p2] = SMP["m_l"]^2; 
SP[l1] = 0; 
SP[l2] = 0; 
SP[k1, k2] = Simplify[(SP[p] - SP[k1] - SP[k2])/2]; 
SP[p, k1] = Simplify[ExpandScalarProduct[SP[k1 + k2, k1]]]; 
SP[p, k2] = Simplify[ExpandScalarProduct[SP[k1 + k2, k2]]]; 
SP[p1, p2] = Simplify[(SP[p] - SP[p1] - SP[p2])/2]; 
SP[p, p1] = Simplify[ExpandScalarProduct[SP[p1 + p2, p1]]]; 
SP[p, p2] = Simplify[ExpandScalarProduct[SP[p1 + p2, p2]]]; 
SP[l1, l2] = Simplify[(SP[p] - SP[l1] - SP[l2])/2]; 
SP[p, l1] = Simplify[ExpandScalarProduct[SP[l1 + l2, l1]]]; 
SP[p, l2] = Simplify[ExpandScalarProduct[SP[l1 + l2, l2]]]; 
```

## Square the amplitudes

We average over the polarizations of the W-boson, hence the additional factor 1/3

```mathematica
ampSquaredLeptonsMassless[0] = 
   Simplify[(DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       DiracSimplify[FermionSpinSum[ampLeptonsMassless[0]*
             ComplexConjugate[ampLeptonsMassless[0]]]]]]
```

![0qmtoh4qwxaud](img/0qmtoh4qwxaud.svg)

```mathematica
ampSquaredLeptonsMassive[0] = 
   Simplify[(DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       DiracSimplify[FermionSpinSum[ampLeptonsMassive[0]*
             ComplexConjugate[ampLeptonsMassive[0]]]]]]
```

![00kauboez11w8](img/00kauboez11w8.svg)

```mathematica
ampSquaredUpQuarks[0] = Simplify[
     (DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       SUNSimplify[DiracSimplify[FermionSpinSum[
             ampUpQuarks[0]*ComplexConjugate[ampUpQuarks[0]]]]]]]
```

![1mrozpl0hw0rm](img/1mrozpl0hw0rm.svg)

```mathematica
ampSquaredDownQuarks[0] = 
   Simplify[(DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       SUNSimplify[DiracSimplify[FermionSpinSum[
             ampDownQuarks[0]*ComplexConjugate[ampDownQuarks[0]]]]]]]
```

![1ws8w1q013wps](img/1ws8w1q013wps.svg)

## Total decay rates

```mathematica
phaseSpacePrefactor[m1_, m2_, M_] := (1/(16*Pi*M))*
       Sqrt[1 - (m1 + m2)^2/M^2]*Sqrt[1 - (m1 - m2)^2/M^2]; 
```

```mathematica
totalDecayRateLeptonsMassless = 
   Simplify[phaseSpacePrefactor[0, 0, SMP["m_Z"]]*
       ampSquaredLeptonsMassless[0]]
```

![1hr8b1lt48lxj](img/1hr8b1lt48lxj.svg)

```mathematica
totalDecayRateLeptonsMassive = 
   Simplify[phaseSpacePrefactor[SMP["m_l"], SMP["m_l"], 
         SMP["m_Z"]]*ampSquaredLeptonsMassive[0]]
```

![0ch2t5l6eeaip](img/0ch2t5l6eeaip.svg)

```mathematica
totalDecayRateUpQuarks = Simplify[
     phaseSpacePrefactor[SMP["m_q"], SMP["m_q"], SMP["m_Z"]]*
       ampSquaredUpQuarks[0]]
```

![0pqoyynahyt6e](img/0pqoyynahyt6e.svg)

```mathematica
totalDecayRateDownQuarks = 
   Simplify[phaseSpacePrefactor[SMP["m_q"], SMP["m_q"], 
         SMP["m_Z"]]*ampSquaredDownQuarks[0]]
```

![1oswkzpoh2bay](img/1oswkzpoh2bay.svg)

## Check the final results

```mathematica
tmp = (SMP["G_F"]*Sqrt[1 - (4*mf^2)/SMP["m_Z"]^2]*SMP["m_Z"]^3*
          (cv^2*(1 + 2*(mf^2/SMP["m_Z"]^2)) + 
             ca^2*(1 - 4*(mf^2/SMP["m_Z"]^2))))/(6*Pi*Sqrt[2]); 
knownResults = {tmp /. {cv | ca -> 1/2, mf -> 0}, 
       tmp /. {cv -> -2^(-1) + 2*SMP["sin_W"]^2, ca -> -2^(-1), 
           mf -> SMP["m_l"]}, CA*tmp /. 
         {cv -> 1/2 - (4/3)*SMP["sin_W"]^2, ca -> 1/2, 
           mf -> SMP["m_q"]}, CA*tmp /. 
         {cv -> -2^(-1) + (2/3)*SMP["sin_W"]^2, ca -> -2^(-1), 
           mf -> SMP["m_q"]}}; 
FCCompareResults[{totalDecayRateLeptonsMassless, 
       totalDecayRateLeptonsMassive, totalDecayRateUpQuarks, 
       totalDecayRateDownQuarks}, knownResults, 
     Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
    Physics, Chapter 5.2:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1k236oo3ls8mk](img/1k236oo3ls8mk.svg)

![1jmwt8nnpjfj9](img/1jmwt8nnpjfj9.svg)