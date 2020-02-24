---
title: Compton scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Ga -> El Ga, QED, matrix element squared, tree"; 
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
       {F[2, {1}], V[1]} -> {F[2, {1}], V[1]}, 
       InsertionLevel -> {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0ka52lbu6hjgo](img/0ka52lbu6hjgo.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, k1}, OutgoingMomenta -> {p2, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1, k2}, List -> False, 
     SMP -> True, Contract -> True]
```

![178s23x7c1i48](img/178s23x7c1i48.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, k1, -p2, -k2, SMP["m_e"], 0, 
     SMP["m_e"], 0]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_e"]^2}] & )[
       DiracSimplify[(FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           (DoPolarizationSums[#1, k2, 0] & )[
             (DoPolarizationSums[#1, k1, 0] & )[
               FeynAmpDenominatorExplicit[amp[0]*ComplexConjugate[
                     amp[0]]]]]]]]]
```

![0t3qenc84j6hq](img/0t3qenc84j6hq.svg)

## Check the final results

```mathematica
knownResult = 2*SMP["e"]^4*(SP[p1, k2]/SP[p1, k1] + 
          SP[p1, k1]/SP[p1, k2] + 2*SMP["m_e"]^2*
            (1/SP[p1, k1] - 1/SP[p1, k2]) + SMP["m_e"]^4*
            (1/SP[p1, k1] - 1/SP[p1, k2])^2); 
FCCompareResults[ampSquared[0], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eq 5.87:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1p5xo4u3xqaze](img/1p5xo4u3xqaze.svg)

![19i94qn3rlt9m](img/19i94qn3rlt9m.svg)