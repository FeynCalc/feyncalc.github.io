---
title: Quark-antiquark pair annihilation into photons
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Q Qbar -> Ga Ga, QCD, matrix element squared, tree"; 
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
       {F[3, {1}], -F[3, {1}]} -> {V[1], V[1]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![01ttj2ymecosd](img/01ttj2ymecosd.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1, k2}, List -> False, 
     SMP -> True, Contract -> True, DropSumOver -> True, 
     Prefactor -> (9/4)*SMP["e_Q"]^2]
```

![0r2tyulwiehgn](img/0r2tyulwiehgn.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_u"], 
     SMP["m_u"], 0, 0]; 
```

## Square the amplitude

We average over the spins and the colors of the quarks, hence the additional factor 1/N^2*1/2^2. 
Since the final state particles are indistinguishable, we add an extra 1/2

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2}] & )[
       (DoPolarizationSums[#1, k2, k1] & )[
         (DoPolarizationSums[#1, k1, k2] & )[DiracSimplify[
             (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
               (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                        False] & )[FeynAmpDenominatorExplicit[
                   (1/2)*(1/SUNN^2)*(amp[0]*ComplexConjugate[
                          amp[0]])]]]]]]]]
```

![1d0lrwg1q71l5](img/1d0lrwg1q71l5.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_u"] -> 0} & )[ampSquared[0]]]
```

![06by8xjbeiwza](img/06by8xjbeiwza.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![103bk42mvg91r](img/103bk42mvg91r.svg)

## Check the final results

```mathematica
knownResults = {((t^2 + u^2)*SMP["e"]^4*SMP["e_Q"]^4)/(3*t*u)}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0]}, {knownResults}, 
   Text -> {"\tCompare to CalcHEP:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![10a8q24v31xe6](img/10a8q24v31xe6.svg)

![18yxvd1lx34ba](img/18yxvd1lx34ba.svg)

![0q1vk9zqdg60e](img/0q1vk9zqdg60e.svg)