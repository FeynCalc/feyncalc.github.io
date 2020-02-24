---
title: Quark-quark scattering (different flavors)
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Qi Qj -> Qi Qj, QCD, matrix element squared, tree"; 
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
       {F[3, {1}], F[3, {2}]} -> {F[3, {1}], F[3, {2}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD", 
       ExcludeParticles -> {S[_], V[1 | 2]}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0deotdkoetijr](img/0deotdkoetijr.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     DropSumOver -> True]
```

![18l762vtopcr4](img/18l762vtopcr4.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_u"], 
     SMP["m_c"], SMP["m_u"], SMP["m_c"]]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2 + 
                2*SMP["m_c"]^2}] & )[DiracSimplify[
         (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                    False] & )[FeynAmpDenominatorExplicit[
               (1/SUNN^2)*(amp[0]*ComplexConjugate[amp[0]])]]]]]]
```

![001bspb65vr7l](img/001bspb65vr7l.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_u" | "m_c"] -> 0} & )[ampSquared[0]]]
```

![0krk1160cssoy](img/0krk1160cssoy.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![07ajwlbjmsjms](img/07ajwlbjmsjms.svg)

## Check the final results

```mathematica
knownResults = {(4/9)*SMP["g_s"]^4*((s^2 + u^2)/t^2)}; 
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

![0g4dnl5ch0b16](img/0g4dnl5ch0b16.svg)