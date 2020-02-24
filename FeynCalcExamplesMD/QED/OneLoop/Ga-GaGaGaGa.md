---
title: 5-photon interaction in QED
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga -> Ga Ga Ga Ga, QED, amplitude, 1-loop"; 
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
MakeBoxes[k4, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(4\)]\)"; 
MakeBoxes[k5, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(5\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 4], 
       {V[1]} -> {V[1], V[1], V[1], V[1]}, InsertionLevel -> 
         {Particles}, ExcludeParticles -> {S[_], V[_], U[_], 
           F[3 | 4], F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {4, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![103dq796unirj](img/103dq796unirj.svg)

![03kz7bihujki7](img/03kz7bihujki7.svg)

![1kg6o3waklmj4](img/1kg6o3waklmj4.svg)

![0eysyl6gqhr0z](img/0eysyl6gqhr0z.svg)

![07qx8n65wnmyf](img/07qx8n65wnmyf.svg)

![0xjodistknmni](img/0xjodistknmni.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, PreFactor -> 1, 
       Truncated -> True], IncomingMomenta -> {k1}, 
     OutgoingMomenta -> {k2, k3, k4, k5}, LoopMomenta -> {q}, 
     LorentzIndexNames -> {mu, nu, rho}, UndoChiralSplittings -> 
       True, ChangeDimension -> D, List -> True, SMP -> True]
```

![0g4lr5jjxv9ea](img/0g4lr5jjxv9ea.svg)

## Calculate the amplitude

We obtain 24 diagrams. The sum vanishes because the contribution of each odd diagram is exactly cancelled by the contribution of the next even diagram, i.e. A1+A2=0, A3+A4=0 and so on

```mathematica
amp[1] = FCTraceFactor[amp[0]]; 
```

```mathematica
amp[2] = Total /@ Partition[amp[1], 2, 2, 1, {}]
```

![1nddb1bbmatog](img/1nddb1bbmatog.svg)

```mathematica
amp[3] = Total[amp[2]]
```

![1q5gsey83nst9](img/1q5gsey83nst9.svg)

## Check the final results

```mathematica
FCCompareResults[amp[3], 0, 
     Text -> 
       {"\tVerify Furry's theorem for 5-photons at 1-loop:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1up54m1dkdg5j](img/1up54m1dkdg5j.svg)

![0zcke7ng2tosq](img/0zcke7ng2tosq.svg)