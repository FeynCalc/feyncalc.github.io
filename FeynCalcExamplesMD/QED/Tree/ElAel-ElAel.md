---
title: Bhabha scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Ael -> El Ael, QED, matrix element squared, tree"; 
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
       {F[2, {1}], -F[2, {1}]} -> {F[2, {1}], -F[2, {1}]}, 
       InsertionLevel -> {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0zhajstebrd2a](img/0zhajstebrd2a.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True]
```

![008979jvg13h7](img/008979jvg13h7.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"], 
     SMP["m_e"], SMP["m_e"], SMP["m_e"]]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
         FeynAmpDenominatorExplicit[amp[0]*ComplexConjugate[
               amp[0]]]]]]
```

![1ktfwgqhypga9](img/1ktfwgqhypga9.svg)

```mathematica
ampSquaredMassless[0] = Simplify[(#1 /. {SMP["m_e"] -> 0} & )[
       ampSquared[0]]]
```

![0h24y9tr85g12](img/0h24y9tr85g12.svg)

## Check the final results

```mathematica
knownResult = 2*SMP["e"]^4*((s^2 + u^2)/t^2) + 
       4*SMP["e"]^4*(u^2/(s*t)) + 2*SMP["e"]^4*((t^2 + u^2)/s^2); 
FCCompareResults[ampSquaredMassless[0], knownResult, 
     Text -> {"\tCheck the final result:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1x113npaxkqx6](img/1x113npaxkqx6.svg)

![0ay1cw1r60e7w](img/0ay1cw1r60e7w.svg)