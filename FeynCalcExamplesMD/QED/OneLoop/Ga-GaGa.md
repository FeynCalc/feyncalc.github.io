---
title: 3-photon interaction in QED
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga -> Ga Ga, QED, amplitude, 1-loop"; 
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
MakeBoxes[nu, TraditionalForm] := "\[Nu]"; 
MakeBoxes[rho, TraditionalForm] := "\[Rho]"; 
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
MakeBoxes[k3, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(3\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 2], 
       {V[1]} -> {V[1], V[1]}, InsertionLevel -> {Particles}, 
       ExcludeParticles -> {S[_], V[_], U[_], F[3 | 4], 
           F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0ztokgr3vlix9](img/0ztokgr3vlix9.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, PreFactor -> 1, 
       Truncated -> True], IncomingMomenta -> {k1}, 
     OutgoingMomenta -> {k2, k3}, LoopMomenta -> {q}, 
     LorentzIndexNames -> {mu, nu, rho}, UndoChiralSplittings -> 
       True, ChangeDimension -> D, List -> False, SMP -> True]
```

![0o8qf3ipebe8j](img/0o8qf3ipebe8j.svg)

## Calculate the amplitude

We obtain two triangle diagrams. The sum vanishes because the contribution of the first diagram cancels the contribution of the second diagram.

```mathematica
amp[1] = FCTraceFactor[amp[0]]
```

![1wob216ap9a8j](img/1wob216ap9a8j.svg)

## Check the final results

```mathematica
FCCompareResults[amp[1], 0, 
     Text -> 
       {"\tVerify Furry's theorem for 3-photons at 1-loop:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![03fm9jbm8icv4](img/03fm9jbm8icv4.svg)

![180jviqcnmhao](img/180jviqcnmhao.svg)