---
title: Quark-antiquark pair production from gluon-gluon annihilation
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Gl Gl -> Q Qbar, QCD, matrix element squared, tree"; 
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
       {V[5], V[5]} -> {F[3, {1}], -F[3, {1}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1778xkvu7vpkv](img/1778xkvu7vpkv.svg)

![10x4iuglrn1v6](img/10x4iuglrn1v6.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {p1, p2}, List -> False, 
     SMP -> True, Contract -> True, DropSumOver -> True]
```

![088p6n6jrec2i](img/088p6n6jrec2i.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, 0, 0, SMP["m_u"]^2, 
     SMP["m_u"]^2]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2}] & )[
       (DoPolarizationSums[#1, p2, p1, ExtraFactor -> 1/2] & )[
         (DoPolarizationSums[#1, p1, p2, ExtraFactor -> 1/2] & )[
           DiracSimplify[FermionSpinSum[
               (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                        False] & )[FeynAmpDenominatorExplicit[
                   (1/(SUNN^2 - 1)^2)*(amp[0]*ComplexConjugate[
                          amp[0]])]]]]]]]]
```

![1hh20sstbaizi](img/1hh20sstbaizi.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_u"] -> 0} & )[ampSquared[0]]]
```

![0gxyg78eg88rx](img/0gxyg78eg88rx.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![098md2c299wei](img/098md2c299wei.svg)

## Check the final results

```mathematica
knownResults = {(1/6)*SMP["g_s"]^4*((t^2 + u^2)/(t*u)) - 
         (3/8)*SMP["g_s"]^4*((t^2 + u^2)/s^2)}; 
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

![0dnlif79tsgg2](img/0dnlif79tsgg2.svg)