---
title: Phi Phi scattering at 1-loop
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom Phi^4 model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/Phi4/GenerateModelPhi4.m before running it for the first time.

```mathematica
description = 
     "Phi Phi -> Phi Phi, Phi^4, asymptotic limit, 1-loop"; 
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

## Configure some options

```mathematica
FAPatch[PatchModelsOnly -> True]; 
```

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
diags = InsertFields[CreateTopologies[1, 2 -> 2, 
         ExcludeTopologies -> {WFCorrections}], 
       {S[1], S[1]} -> {S[1], S[1]}, InsertionLevel -> {Classes}, 
       Model -> FileNameJoin[{"Phi4", "Phi4"}]]; 
Paint[diags, ColumnsXRows -> {3, 1}, Numbering -> None, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0oeasj2l0zgbe](img/0oeasj2l0zgbe.svg)

```mathematica
diagsCT = InsertFields[CreateCTTopologies[1, 2 -> 2, 
         ExcludeTopologies -> {WFCorrectionCTs}], 
       {S[1], S[1]} -> {S[1], S[1]}, InsertionLevel -> {Classes}, 
       Model -> FileNameJoin[{"Phi4", "Phi4"}]]; 
Paint[diagsCT, ColumnsXRows -> {1, 1}, Numbering -> None, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![0me6m948a5kjm](img/0me6m948a5kjm.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, PreFactor -> 1], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     LoopMomenta -> {q}, ChangeDimension -> D, List -> False, 
     FinalSubstitutions -> {Mphi -> m}]
```

![1hkcvfg4hs5h1](img/1hkcvfg4hs5h1.svg)

```mathematica
ampCT[0] = FCFAConvert[CreateFeynAmp[diagsCT, PreFactor -> 1], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     LoopMomenta -> {q}, ChangeDimension -> D, List -> False, 
     FinalSubstitutions -> {Mphi -> m, 
         Zg -> 1 + SMP["d_g^MSbar"]}]
```

![0wdhj6l4ikgu7](img/0wdhj6l4ikgu7.svg)

## Fix the kinematics

For simplicity, let us consider the massless case

```mathematica
FCClearScalarProducts[]
SetMandelstam[s, t, u, p1, p2, -k1, -k2, 0, 0, 0, 0]; 
```

## Calculate the amplitude

```mathematica
amp[1] = (ToPaVe[#1, q] & )[(#1 /. m -> 0 & )[amp[0]]]
```

![0ouydg2g7dyak](img/0ouydg2g7dyak.svg)

The explicit value of the integral can be obtained from Package-X via the FeynHelpers add-on.

```mathematica
loopInt = {B0[s_, 0, 0] :> 
         -(-2 + Log[4*Pi] - Log[(-4*Pi*ScaleMu^2)/s])/(16*Pi^4) + 
           SMP["Delta"]/(16*Pi^4)}; 
```

```mathematica
amp[2] = Simplify[amp[1] /. loopInt]
```

![0vnuiqz9nyc93](img/0vnuiqz9nyc93.svg)

```mathematica
ampFull[0] = Expand[amp[2] + ampCT[0] /. SMP["d_g^MSbar"] -> 
         (3*g^2*SMP["Delta"])/(32*Pi^2)]
```

![1ivc1qkl2ezer](img/1ivc1qkl2ezer.svg)

```mathematica
FCCompareResults[FreeQ[ampFull, SMP["Delta"]], True, 
     Text -> 
       {"\tThe UV divergence is cancelled by the counter-term:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
```

![1580xnzlcztjx](img/1580xnzlcztjx.svg)

Now let us look at the asymptotic limit where  s goes to infinity and t is fixed

```mathematica
ampFullAsy[0] = Normal[Series[ampFull[0] /. u -> -s - t, 
       {s, Infinity, 0}]]
```

![0q3ablo9zyty4](img/0q3ablo9zyty4.svg)

The leading order behavior is governed by the log of s

```mathematica
ampFullAsy[1] = (SelectNotFree2[#1, s] & )[
     PowerExpand[ampFullAsy[0]]]
```

![04303zvz8v6en](img/04303zvz8v6en.svg)

## Check the final results

```mathematica
knownResult = ((-I/16)*g^2*Log[s])/Pi^2; 
FCCompareResults[ampFullAsy[1], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Ex 10.4:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1dowyj2d5odnw](img/1dowyj2d5odnw.svg)

![13rk8r467txo3](img/13rk8r467txo3.svg)