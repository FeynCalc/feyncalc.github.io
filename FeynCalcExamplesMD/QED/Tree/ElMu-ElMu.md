---
title: Electron-muon scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Mu -> El Mu, QED, matrix element squared, tree"; 
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
       {F[2, {1}], F[2, {2}]} -> {F[2, {1}], F[2, {2}]}, 
       InsertionLevel -> {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![18f208060zjpr](img/18f208060zjpr.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True]
```

![05je5u69suo0h](img/05je5u69suo0h.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"], 
     SMP["m_mu"], SMP["m_e"], SMP["m_mu"]]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
         FeynAmpDenominatorExplicit[amp[0]*ComplexConjugate[
               amp[0]]]]]]
```

![1g95bi5qkafkv](img/1g95bi5qkafkv.svg)

```mathematica
ampSquaredMassless1[0] = Simplify[(#1 /. {SMP["m_e"] -> 0} & )[
       ampSquared[0]]]
```

![1dsz8v7hrul5p](img/1dsz8v7hrul5p.svg)

```mathematica
ampSquaredMassless2[0] = Simplify[
     (#1 /. {SMP["m_e"] -> 0, SMP["m_mu"] -> 0} & )[
       ampSquared[0]]]
```

![0ue52ozqzy1dq](img/0ue52ozqzy1dq.svg)

## Check the final results

```mathematica
knownResults = {(#1 /. SMP["m_e"] -> 0 & )[
         ExpandScalarProduct[(8*SMP["e"]^4*(SP[p1, k2]*SP[p2, k1] + 
                   SP[p1, p2]*SP[k1, k2] - 
          SMP["m_mu"]^2*SP[p1, k1]))/
             SP[k1 - p1]^2]], (8*(SMP["e"]^4/t^2))*
         ((s/2)^2 + (u/2)^2)}; 
FCCompareResults[{ampSquaredMassless1[0], 
     ampSquaredMassless2[0]}, knownResults, 
   Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
   to QFT, Eqs 5.61 and 5.71:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1j4583ob4y69j](img/1j4583ob4y69j.svg)

![18yxvd1lx34ba](img/18yxvd1lx34ba.svg)

![0qsq0h70rouzq](img/0qsq0h70rouzq.svg)