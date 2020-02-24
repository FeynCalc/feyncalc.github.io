---
title: QCD ghost self-energy
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Gh -> Gh, QCD, only UV divergences, 1-loop"; 
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

We keep scaleless B0 functions, since otherwise the UV part would not come out right.

```mathematica
$KeepLogDivergentScalelessIntegrals = True; 
```

## Generate Feynman diagrams

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 1, 
         ExcludeTopologies -> {Tadpoles}], {U[5]} -> {U[5]}, 
       InsertionLevel -> {Particles}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![1l360erikgn5p](img/1l360erikgn5p.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
       GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p}, LoopMomenta -> {q}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, DropSumOver -> True, 
     Contract -> True]
```

![1d0yu9owfagax](img/1d0yu9owfagax.svg)

## Calculate the amplitude

```mathematica
amp[1] = (TID[#1, q, ToPaVe -> True] & )[SUNSimplify[amp[0]]]
```

![10kq8sh3roxbg](img/10kq8sh3roxbg.svg)

The UV divergence of the amplitude can be obtained via PaVeUVPart.
Here we also need to reintroduce the implicit 1/(2Pi)^D prefactor.
Hint: If you need the full result for the amplitude, use PaXEvaluate from FeynHelpers.

```mathematica
ampDiv[0] = Simplify[(SelectNotFree2[#1, Epsilon] & )[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
             PaVeUVPart[amp[1], Prefactor -> 1/(2*Pi)^D]]]]]]
```

![1hvj9cu7qlpiv](img/1hvj9cu7qlpiv.svg)

The self-energy amplitude is usually defined as  (p^2 delta^ab  Pi(p^2)

```mathematica
pi[0] = Cancel[FCI[ampDiv[0]/(I*SUNDelta[SUNIndex[Glu1], 
              SUNIndex[Glu2]]*SPD[p, p])]]
```

![1f5zm85mwwfzg](img/1f5zm85mwwfzg.svg)

## Check the final results

```mathematica
knownResult = (-SMP["g_s"]^2/(4*Pi)^2)*CA*((3 - GaugeXi[g])/4)*
       (1/Epsilon); 
FCCompareResults[pi[0], knownResult, 
     Text -> 
       {"\tCompare to Muta, Foundations of QCD, Eq. 2.5.136:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0p3bsv1nbiiw7](img/0p3bsv1nbiiw7.svg)

![1wiyh3vtaplaq](img/1wiyh3vtaplaq.svg)