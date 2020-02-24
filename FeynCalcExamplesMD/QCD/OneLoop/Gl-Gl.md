---
title: QCD vacuum polarization
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Gl -> Gl, QCD, only UV divergences, 1-loop"; 
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
diags = InsertFields[CreateTopologies[1, 1 -> 1, 
         ExcludeTopologies -> {Tadpoles}], {V[5]} -> {V[5]}, 
       InsertionLevel -> {Particles}, Model -> "SMQCD", 
       ExcludeParticles -> {S[_], V[2 | 3], F[4], F[3, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {2, 2}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 512}]; 
```

![1rq4lch23fya4](img/1rq4lch23fya4.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
       GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p}, LoopMomenta -> {q}, 
     LorentzIndexNames -> {mu, nu}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> True, SMP -> True, 
     DropSumOver -> True, Contract -> True, 
     FinalSubstitutions -> {SMP["m_u"] -> SMP["m_q"]}]
```

![0chkdse8cah68](img/0chkdse8cah68.svg)

## Calculate the amplitude

### The gluon tadpole 

This contribution is zero in dimensional regularization, because the loop integrals have no scale (and they are not log divergent)

```mathematica
amp1[0] = TID[amp[0][[1]], q, ToPaVe -> True]
```

![1rtj0hgv0ebfl](img/1rtj0hgv0ebfl.svg)

```mathematica
FCCompareResults[amp1[0], 0, 
     Text -> {"\tThe gluon tadpole vanishes:", "CORRECT.", 
         "WRONG!"}, Interrupt -> {Hold[Quit[1]], Automatic}]; 
```

![0wfajb5pcz0rd](img/0wfajb5pcz0rd.svg)

### The quark loop

```mathematica
amp2[0] = (TID[#1, q, ToPaVe -> True] & )[
     SUNSimplify[amp[0][[2]]]]
```

![0qeqvn3gx4eqy](img/0qeqvn3gx4eqy.svg)

The contribution of the quark loop alone is  gauge invariant.

```mathematica
tmp = Factor[Contract[FVD[p, mu]*FVD[p, nu]*amp2[0]]]
FCCompareResults[tmp, 0, Text -> 
       {"\tThe quark loop contribution is gauge invariant:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
```

![1wob216ap9a8j](img/1wob216ap9a8j.svg)

![14g332x5t9c2h](img/14g332x5t9c2h.svg)

### The ghost loop

```mathematica
amp3[0] = (TID[#1, q, ToPaVe -> True] & )[
     SUNSimplify[amp[0][[3]]]]
```

![0z3o0j6hw3s1i](img/0z3o0j6hw3s1i.svg)

The contribution of the gluon loop alone is not gauge invariant.

```mathematica
tmp1 = Factor[Contract[FVD[p, mu]*FVD[p, nu]*amp3[0]]]
```

![0vgi4seb6e6sx](img/0vgi4seb6e6sx.svg)

### The gluon loop

```mathematica
amp4[0] = (TID[#1, q, ToPaVe -> True] & )[
     SUNSimplify[amp[0][[4]]]]
```

![1nuoiwahcj5fy](img/1nuoiwahcj5fy.svg)

The contribution of the gluon loop alone is not gauge invariant. Notice, however, that the sum
of the ghost and gluon contributions is gauge invariant!

```mathematica
tmp2 = Factor[Contract[FVD[p, mu]*FVD[p, nu]*amp4[0]]]
```

![1gom845ahm6ba](img/1gom845ahm6ba.svg)

```mathematica
FCCompareResults[tmp1 + tmp2, 0, 
     Text -> {"\tThe sum of the ghost and gluon loop \
    contributions is gauge invariant:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
```

![1osdhpwfexkrc](img/1osdhpwfexkrc.svg)

### Putting everything together

When adding all the contributions together, we multiply the quark contribution by N_f to account for the 6 quark flavors that actually run in that loop. We ignore the fact that different flavors have different masses, since the divergent piece of the gluon self-energy will not depend on the quark mass.

```mathematica
amp[1] = Nf*amp2[0] + amp3[0] + amp4[0]
```

![0ymefksb8iujo](img/0ymefksb8iujo.svg)

The UV divergence of the amplitude can be obtained via PaVeUVPart.
Here we also need to reintroduce the implicit 1/(2Pi)^D prefactor.
Hint: If you need the full result for the amplitude, use PaXEvaluate from FeynHelpers.

```mathematica
ampDiv[0] = Simplify[(SelectNotFree2[#1, Epsilon] & )[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
             PaVeUVPart[amp[1], Prefactor -> 1/(2*Pi)^D]]]]]]
```

![1mtc6a8vcgz5a](img/1mtc6a8vcgz5a.svg)

The self-energy amplitude is usually defined as  (p^2 g^{mu nu} - p^mu p^nu) i Pi(p^2)

```mathematica
pi[0] = Cancel[FCI[ampDiv[0]/(I*SUNDelta[SUNIndex[Glu1], 
              SUNIndex[Glu2]]*(SPD[p, p]*MTD[mu, nu] - 
               FVD[p, mu]*FVD[p, nu]))]]
```

![1a60fiyvs04uc](img/1a60fiyvs04uc.svg)

## Check the final results

```mathematica
knownResult = (-(SMP["g_s"]^2/(4*Pi)^2))*
       ((4/3)*(1/2)*Nf - (1/2)*CA*(13/3 - GaugeXi[g]))*
       (1/Epsilon); 
FCCompareResults[pi[0], knownResult, 
     Text -> {"\tCompare to Muta, Foundations of QCD, Eqs \
    2.5.131-2.5.132:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0467wfifnpi1g](img/0467wfifnpi1g.svg)

![0x0w2veyqf45b](img/0x0w2veyqf45b.svg)