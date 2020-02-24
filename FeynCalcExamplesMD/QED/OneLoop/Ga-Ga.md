---
title: QED vacuum polarization
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga -> Ga, QED, only UV divergences, 1-loop"; 
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

Nicer typesetting

```mathematica
MakeBoxes[mu, TraditionalForm] := "\[Mu]"; 
MakeBoxes[nu, TraditionalForm] := "\[Nu]"; 
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 1], 
       {V[1]} -> {V[1]}, InsertionLevel -> {Particles}, 
       ExcludeParticles -> {S[_], V[2 | 3], (S | U)[_], F[3 | 4], 
           F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![15cq7m7ff5g3b](img/15cq7m7ff5g3b.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
       PreFactor -> 1], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p}, LoopMomenta -> {q}, 
     LorentzIndexNames -> {mu, nu}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> False, SMP -> True, 
     Contract -> True]
```

![12i7qyzie9owg](img/12i7qyzie9owg.svg)

## Calculate the amplitude

```mathematica
amp[1] = TID[amp[0], q, ToPaVe -> True]
```

![0orlwescia6xj](img/0orlwescia6xj.svg)

Check the gauge invariance

```mathematica
tmp = Factor[Contract[FVD[p, mu]*FVD[p, nu]*amp[1]]]
FCCompareResults[tmp, 0, Text -> 
       {"\tThe photon self-energy is gauge invariant:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
```

![1fex7f5pysep6](img/1fex7f5pysep6.svg)

![0gtx1traag30n](img/0gtx1traag30n.svg)

The UV divergence of the amplitude can be obtained via PaVeUVPart.
Here we also need to reintroduce the implicit 1/(2Pi)^D prefactor.
Hint: If you need the full result for the amplitude, use PaXEvaluate from FeynHelpers.

```mathematica
ampDiv[0] = Simplify[(SelectNotFree2[#1, Epsilon] & )[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
             PaVeUVPart[amp[1], Prefactor -> 1/(2*Pi)^D]]]]]]
```

![1vps0ibvsmc2q](img/1vps0ibvsmc2q.svg)

The self-energy amplitude is usually defined as  (p^2 g^{mu nu} - p^mu p^nu) i Pi(p^2)

```mathematica
pi[0] = Cancel[FCI[ampDiv[0]/(I*(SPD[p, p]*MTD[mu, nu] - 
               FVD[p, mu]*FVD[p, nu]))]]
```

![1eevwrbzh3cs5](img/1eevwrbzh3cs5.svg)

## Check the final results

Keep in mind that Peskin and Schroeder use D = 4-Epsilon,
while we did the calculation with D = 4-2Epsilon.

```mathematica
knownResult = (#1 /. 1/Epsilon -> 1/(2*Epsilon) & )[
       (Integrate[#1, {x, 0, 1}] & )[
         (SelectNotFree2[#1, Epsilon] & )[
           Normal[(Series[#1, {Epsilon, 0, 0}] & )[
               (FCReplaceD[#1, D -> 4 - Epsilon] & )[
                 (-SMP["e"]^2/(4*Pi)^(D/2))*(Gamma[2 - D/2]/
                      (SMP["m_e"]^2 - x*(1 - x)*SPD[p, p])^(2 - 
              D/2))*
                   (8*x*(1 - x))]]]]]]; 
FCCompareResults[pi[0], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eq 10.44:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0vrcllxu5kif3](img/0vrcllxu5kif3.svg)

![1rjdrtz6trhck](img/1rjdrtz6trhck.svg)