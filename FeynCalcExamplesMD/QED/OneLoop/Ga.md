---
title: Photon tadpole in QED
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga, QED, amplitude, 1-loop"; 
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
MakeBoxes[mu, TraditionalForm] := "\[Mu]"; 
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 0], {V[1]} -> {}, 
       InsertionLevel -> {Particles}, ExcludeParticles -> 
         {S[_], V[_], U[_], F[3 | 4], F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![1vbor99cu217z](img/1vbor99cu217z.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, PreFactor -> 1, 
       Truncated -> True], IncomingMomenta -> {k}, 
     LorentzIndexNames -> {mu}, LoopMomenta -> {q}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True]
```

![08whtdg4qnt79](img/08whtdg4qnt79.svg)

## Calculate the amplitude

Having performed the Dirac algebra we clearly see that this diagram must 
vanish because the loop integral is antisymmetric under q^mu -> - q^mu.

```mathematica
amp[1] = DiracSimplify[amp[0]]
```

![0plj14iwtfbdy](img/0plj14iwtfbdy.svg)

TID can recognize this and we obtain zero

```mathematica
amp[2] = TID[amp[1], q]
```

![11dm8hxkf5oud](img/11dm8hxkf5oud.svg)

## Check the final results

```mathematica
FCCompareResults[amp[2], 0, 
     Text -> {"\tVerify Furry's theorem for 1-photon at 1-loop:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![11hzb2hvyoz4g](img/11hzb2hvyoz4g.svg)

![0fidso3wnip61](img/0fidso3wnip61.svg)