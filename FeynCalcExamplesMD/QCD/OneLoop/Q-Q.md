---
title: QCD quark self-energy
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Q -> Q, QCD, only UV divergences, 1-loop"; 
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
         ExcludeTopologies -> Tadpoles], {F[3, {1}]} -> 
         {F[3, {1}]}, InsertionLevel -> {Particles}, 
       Model -> "SMQCD", ExcludeParticles -> 
         {S[_], V[1 | 2 | 3]}]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![1628uyfltrg4i](img/1628uyfltrg4i.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit. We keep the full gauge dependence.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
       PreFactor -> 1, GaugeRules -> {}], IncomingMomenta -> {p}, 
     OutgoingMomenta -> {p}, LoopMomenta -> {q}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, DropSumOver -> True, 
     Contract -> True, FinalSubstitutions -> 
       {SMP["m_u"] -> SMP["m_q"]}]
```

![17a8m21hs2t9k](img/17a8m21hs2t9k.svg)

## Calculate the amplitude

```mathematica
amp[1] = (TID[#1, q, ToPaVe -> True] & )[SUNSimplify[amp[0]]]
```

![0oy9stbxh3rtq](img/0oy9stbxh3rtq.svg)

The UV divergence of the amplitude can be obtained via PaVeUVPart.
Here we also need to reintroduce the implicit 1/(2Pi)^D prefactor.
Hint: If you need the full result for the amplitude, use PaXEvaluate from FeynHelpers.

```mathematica
ampDiv[0] = Simplify[(SelectNotFree2[#1, Epsilon] & )[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
             PaVeUVPart[amp[1], Prefactor -> 1/(2*Pi)^D]]]]]]
```

![0i8g3xx3iv6sp](img/0i8g3xx3iv6sp.svg)

The self-energy amplitude is usually defined as -i Sigma(p^2)

```mathematica
sigma[0] = I*ampDiv[0]
```

![1522c1fs4n8ux](img/1522c1fs4n8ux.svg)

```mathematica
sigmaFeynmanGauge[0] = sigma[0] /. GaugeXi[g] -> 1
```

![04g66cksd5gv8](img/04g66cksd5gv8.svg)

## Check the final results

Notice that the result in the book must be multiplied by (-1) due to the way how self-energy is defined there (c.f. Eq. 2.4.4 and Eq. 2.4.6).

```mathematica
knownResult = (-((-SMP["g_s"]^2/(4*Pi)^2)*CF*(3 + GaugeXi[g])*
               (1/Epsilon)*SMP["m_q"] + GSD[p]*(SMP["g_s"]^2/(4*Pi)^2)*
               CF*GaugeXi[g]*(1/Epsilon)))*SDF[Col1, Col2]; 
FCCompareResults[sigma[0], knownResult, 
   Text -> {"\tCompare to Muto, Foundations of QCD, Eq 10.41:", 
       "CORRECT.", "WRONG!"}, Interrupt -> 
     {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0mh4yjjvk09rs](img/0mh4yjjvk09rs.svg)

![12hkzd5e0qf2q](img/12hkzd5e0qf2q.svg)

![13jexjyrrcxlg](img/13jexjyrrcxlg.svg)