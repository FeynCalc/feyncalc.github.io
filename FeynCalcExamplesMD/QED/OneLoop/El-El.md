---
title: QED electron self-energy
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "El -> El, QED, only UV divergences, 1-loop"; 
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
         ExcludeTopologies -> Tadpoles], {F[2, {1}]} -> 
         {F[2, {1}]}, InsertionLevel -> {Particles}, 
       ExcludeParticles -> {S[_], V[2 | 3], (S | U)[_], F[3 | 4], 
           F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![152iik2s27no5](img/152iik2s27no5.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit. We keep the full gauge dependence.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
       PreFactor -> 1, GaugeRules -> {}], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p}, LoopMomenta -> {q}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, Contract -> True]
```

![0ufzezsghtzak](img/0ufzezsghtzak.svg)

## Calculate the amplitude

```mathematica
amp[1] = TID[amp[0], q, ToPaVe -> True]
```

![1o1wikzm55ijb](img/1o1wikzm55ijb.svg)

The UV divergence of the amplitude can be obtained via PaVeUVPart.
Here we also need to reintroduce the implicit 1/(2Pi)^D prefactor.
Hint: If you need the full result for the amplitude, use PaXEvaluate from FeynHelpers.

```mathematica
ampDiv[0] = Simplify[(SelectNotFree2[#1, Epsilon] & )[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
             PaVeUVPart[amp[1], Prefactor -> 1/(2*Pi)^D]]]]]]
```

![04i5d6wtjl3xc](img/04i5d6wtjl3xc.svg)

The self-energy amplitude is usually defined as -i Sigma(p^2)

```mathematica
sigma[0] = I*ampDiv[0]
```

![112pw3pp6bsmy](img/112pw3pp6bsmy.svg)

```mathematica
sigmaFeynmanGauge[0] = sigma[0] /. GaugeXi[A] -> 1
```

![00oohmpenoz4k](img/00oohmpenoz4k.svg)

## Check the final results

Keep in mind that Peskin and Schroeder use D = 4-Epsilon,
while we did the calculation with D = 4-2Epsilon.

```mathematica
knownResult = (#1 /. 1/Epsilon -> 1/(2*Epsilon) & )[
       (Integrate[#1, {x, 0, 1}] & )[
         (SelectNotFree2[#1, Epsilon] & )[
           Normal[(Series[#1, {Epsilon, 0, 0}] & )[
               (FCReplaceD[#1, D -> 4 - Epsilon] & )[
                 (SMP["e"]^2/(4*Pi)^(D/2))*(Gamma[2 - D/2]/
                      ((1 - x)*SMP["m_e"]^2 + x*ScaleMu^2 - 
                           x*(1 - x)*SPD[p, p])^(2 - D/2))*
                   ((4 - Epsilon)*SMP["m_e"] - (2 - Epsilon)*x*
                        GSD[p])]]]]]]; 
FCCompareResults[sigmaFeynmanGauge[0], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eq 10.41:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0o9lwchf1nwe5](img/0o9lwchf1nwe5.svg)

![039oiap9554pq](img/039oiap9554pq.svg)