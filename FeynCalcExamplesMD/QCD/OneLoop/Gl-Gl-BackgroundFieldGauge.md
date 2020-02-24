---
title: Pure Yang-Mills 1-loop gluon self-energy in the background field formalism
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Gl -> Gl, YM+BGF, only UV divergences, 1-loop"; 
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
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 1, 
         ExcludeTopologies -> {Tadpoles}], {V[50, {a}]} -> 
         {V[50, {b}]}, InsertionLevel -> {Classes}, 
       Model -> FileNameJoin[{"QCDBGF", "QCDBGF"}], 
       GenericModel -> FileNameJoin[{"QCDBGF", "QCDBGF"}], 
       ExcludeParticles -> {F[_]}]; 
Paint[diags, ColumnsXRows -> {2, 2}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 512}]; 
```

![1n771ds0xnaji](img/1n771ds0xnaji.svg)

## Obtain corresponding amplitudes

The 1/(2Pi)^D prefactor is implicit.

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, Truncated -> True, 
         GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p}, 
       OutgoingMomenta -> {p}, LoopMomenta -> {l}, 
       LorentzIndexNames -> {mu, nu}, UndoChiralSplittings -> 
         True, ChangeDimension -> D, List -> True, SMP -> True, 
       DropSumOver -> True, FinalSubstitutions -> 
         {SMP["m_u"] -> SMP["m_q"], GaugeXi[V[5, {_}]] :> 
             GaugeXi[G]}]; 
```

```mathematica
amp[1] = DiracSimplify /@ amp[0]; 
```

```mathematica
amp[2] = (SUNSimplify[TID[#1, l, ToPaVe -> True]] & ) /@ 
       amp[1]; 
```

Discard all the finite pieces of the 1-loop amplitude.

```mathematica
ampDiv[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & ) /@ 
     amp[2]
```

![1aezkx289xwy4](img/1aezkx289xwy4.svg)

```mathematica
ampDiv[1] = Simplify[Normal[(Series[#1, {Epsilon, 0, -1}] & )[
         FCReplaceD[ampDiv[0], D -> 4 - 2*Epsilon]]]]
```

![0sp0nxe6wdpyx](img/0sp0nxe6wdpyx.svg)

## Check the final results

```mathematica
knownResult = FCI[{0, 0, I*CA*SMP["g_s"]^2*
           (SUNDelta[a, b]/(4*Pi)^2)*(1/(3*Epsilon))*
           (MTD[mu, nu]*SPD[p] - FVD[p, mu]*FVD[p, nu]), 
         I*CA*SMP["g_s"]^2*(SUNDelta[a, b]/(4*Pi)^2)*
           (10/(3*Epsilon))*(MTD[mu, nu]*SPD[p] - 
              FVD[p, mu]*FVD[p, nu])}]; 
FCCompareResults[ampDiv[1] /. GaugeXi[G] -> 1, knownResult, 
     Text -> {"\tCompare to Abbott, Nucl. Phys. B 185 (1981) \
    189-203, Eqs 5.11-5.12:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0fqig737ysd8b](img/0fqig737ysd8b.svg)

![1h3vj3h73vvq6](img/1h3vj3h73vvq6.svg)