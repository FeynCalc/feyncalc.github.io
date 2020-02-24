---
title: Antiup quark down quark annihilation into an electron and an antielectron-neutrino
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Qutbar Qdt -> Nel Anel, EW, matrix element squared, tree"; 
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

Enable CKM mixing

```mathematica
$CKM = True; 
```

To avoid dealing with Goldstone bosons we do  the computation in the unitary gauge

```mathematica
InitializeModel[{SM, UnitarySM}, GenericModel -> 
       {Lorentz, UnitaryLorentz}]; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {-F[3, {1}], F[4, {1}]} -> {F[2, {1}], -F[1, {1}]}, 
       InsertionLevel -> {Particles}, Model -> {SM, UnitarySM}, 
       GenericModel -> {Lorentz, UnitaryLorentz}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0rrmq3nimp0gy](img/0rrmq3nimp0gy.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, 
       GaugeRules -> {FAGaugeXi[W | Z] -> Infinity}], 
     IncomingMomenta -> {p2, p1}, OutgoingMomenta -> {k1, k2}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True]
```

![0psv421zulksm](img/0psv421zulksm.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]
SetMandelstam[s, t, u, p1, p2, -k1, -k2, 0, 0, 0, 0]; 
```

## Square the amplitude

We average over the spins and the colors of the quarks, hence the additional factor 1/3^2 1/2^2.

```mathematica
ampSquared[0] = (#1 /. SUNN -> 3 & )[
     (SUNSimplify[#1, SUNNToCACF -> False] & )[
       FeynAmpDenominatorExplicit[DiracSimplify[
           (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
             (1/3^2)*(amp[0]*ComplexConjugate[amp[0]])]]]]]
```

![0gxxd1mcxbt1n](img/0gxxd1mcxbt1n.svg)

## Check the final results

```mathematica
knownResults = {(u^2*SMP["e"]^4*SMP["V_ud", -I]*SMP["V_ud", I])/
         (12*(s - SMP["m_W"]^2)^2*SMP["sin_W"]^4)}; 
FCCompareResults[{ampSquared[0]}, knownResults, 
     Text -> {"\tCompare to CompHEP:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1vx8chj5d7mbj](img/1vx8chj5d7mbj.svg)

![0t7ob9roo34ep](img/0t7ob9roo34ep.svg)