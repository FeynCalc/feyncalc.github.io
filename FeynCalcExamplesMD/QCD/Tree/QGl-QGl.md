---
title: Quark-gluon scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Q Gl -> Q Gl, QCD, matrix element squared, tree"; 
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
       {F[3, {1}], V[5]} -> {F[3, {1}], V[5]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 2}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 512}]; 
```

![0ef7sqb0l9jsa](img/0ef7sqb0l9jsa.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, k1}, OutgoingMomenta -> {p2, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1, k2}, List -> False, 
     SMP -> True, Contract -> True, DropSumOver -> True]
```

![1u43st7ex0i4p](img/1u43st7ex0i4p.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, k1, -p2, -k2, SMP["m_u"], 0, 
     SMP["m_u"], 0]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2}] & )[
       (DoPolarizationSums[#1, k2, k1] & )[
         (DoPolarizationSums[#1, k1, k2, ExtraFactor -> 1/2] & )[
           DiracSimplify[(FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
               (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                        False] & )[FeynAmpDenominatorExplicit[
                   (1/(SUNN*(SUNN^2 - 1)))*(amp[0]*ComplexConjugate[
                          amp[0]])]]]]]]]]
```

![110izgrjxw1e4](img/110izgrjxw1e4.svg)

```mathematica
ampSquaredMassless[0] = Simplify[(#1 /. {SMP["m_u"] -> 0} & )[
       ampSquared[0]]]
```

![19z5y4nhyvcq5](img/19z5y4nhyvcq5.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![1q7zqfdbkih9d](img/1q7zqfdbkih9d.svg)

## Check the final results

```mathematica
knownResults = {(-(4/9))*SMP["g_s"]^4*((s^2 + u^2)/(s*u)) + 
         SMP["g_s"]^4*((u^2 + s^2)/t^2)}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0]}, {knownResults}, 
   Text -> {"\tCompare to Ellis, Stirling and Weber, QCD and \
   Collider Physics, Table 7.1:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}, 
   Factoring -> Function[x, Simplify[TrickMandelstam[x, 
           {s, t, u, 0}]]]]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1sjri13b79q9g](img/1sjri13b79q9g.svg)

![18yxvd1lx34ba](img/18yxvd1lx34ba.svg)

![02g9p8jc1q2ap](img/02g9p8jc1q2ap.svg)