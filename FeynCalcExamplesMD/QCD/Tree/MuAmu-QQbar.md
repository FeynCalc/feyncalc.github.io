---
title: Quark-antiquark production in muon-antimuon annihilation
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Mu Amu -> Q Qbar, QCD, total cross section, tree"; 
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
       {F[2, {2}], -F[2, {2}]} -> {F[3, {1}], -F[3, {1}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD", 
       ExcludeParticles -> {S[_], V[2]}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1lodhonf4s1mo](img/1lodhonf4s1mo.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     DropSumOver -> True, Prefactor -> (3/2)*SMP["e_Q"]]
```

![07o6911z9sv0y](img/07o6911z9sv0y.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_mu"]^2, 
     SMP["m_mu"]^2, SMP["m_u"]^2, SMP["m_u"]^2]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2 + 
                2*SMP["m_mu"]^2}] & )[DiracSimplify[
         (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                    False] & )[FeynAmpDenominatorExplicit[
               amp[0]*ComplexConjugate[amp[0]]]]]]]]
```

![1r1lchsjivob1](img/1r1lchsjivob1.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_u" | "m_mu"] -> 0} & )[ampSquared[0]]]
```

![0mlek221ua9tu](img/0mlek221ua9tu.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![11lsfmp61x82y](img/11lsfmp61x82y.svg)

## Total cross-section

```mathematica
integral = Integrate[Simplify[ampSquaredMasslessSUNN3[0]/
             (s/4) /. u -> -s - t], {t, -s, 0}] /. 
     SMP["e"]^4 -> (4*Pi*SMP["alpha_fs"])^2
```

![1njz58dxub28w](img/1njz58dxub28w.svg)

```mathematica
prefac = 2*(Pi/(128*Pi^2*s))
```

![14ddjetrwomxi](img/14ddjetrwomxi.svg)

The total cross-section 

```mathematica
crossSectionTotal = Factor2[PowerExpand[integral*prefac]]
```

![17sm4r5lthaz2](img/17sm4r5lthaz2.svg)

## Check the final results

```mathematica
knownResults = {(6*(t^2 + u^2)*SMP["e"]^4*SMP["e_Q"]^2)/s^2, 
       (4*Pi*SMP["alpha_fs"]^2*SMP["e_Q"]^2)/s}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0], 
     crossSectionTotal}, knownResults, 
   Text -> {"\tCompare to CalcHEP and to Field, Applications of \
   Perturbative QCD, Eq. 2.1.15", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1uxwq9qvk4ajk](img/1uxwq9qvk4ajk.svg)

![1311eokuj92n7](img/1311eokuj92n7.svg)

![11oxdbyn8nacs](img/11oxdbyn8nacs.svg)