---
title: 1-loop phi^3 renormalization in the minimal subtraction schemes
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom phi^3 model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/Phi3/GenerateModelPhi3.m before running it for the first time.

```mathematica
description = "Renormalization, phi^3, MS and MSbar, 1-loop"; 
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

```mathematica
FAPatch[PatchModelsOnly -> True]; 
```

## Generate Feynman diagrams

```mathematica
params = {InsertionLevel -> {Particles}, 
       Model -> FileNameJoin[{"Phi3", "Phi3"}], 
       GenericModel -> FileNameJoin[{"Phi3", "Phi3"}]}; 
top[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> Tadpoles]; 
topCT[i_, j_] := CreateCTTopologies[1, i -> j]; 
topVertex[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> {WFCorrections}]; 
topVertexCT[i_, j_] := CreateCTTopologies[1, i -> j, 
       ExcludeTopologies -> {WFCorrectionCTs}]; 
{diagPhi3SE, diagPhi3SECT} = 
     (InsertFields[#1, {S[1]} -> {S[1]}, Sequence @@ 
              params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagVertex, diagVertexCT} = 
     (InsertFields[#1, {S[1]} -> {S[1], S[1]}, 
            Sequence @@ params] & ) /@ {topVertex[1, 2], 
         topVertexCT[1, 2]}; 
```

```mathematica
diag1[0] = diagPhi3SE[[0]][Sequence @@ diagPhi3SE, 
       Sequence @@ diagPhi3SECT]; 
diag2[0] = diagVertex[[0]][Sequence @@ diagVertex, 
       Sequence @@ diagVertexCT]; 
```

```mathematica
Paint[diag1[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![1m7xanqldq8lh](img/1m7xanqldq8lh.svg)

```mathematica
Paint[diag2[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![0oi1pluqwcw1j](img/0oi1pluqwcw1j.svg)

## Obtain the amplitudes

The 1/(2Pi)^D prefactor is implicit.

Self-energy including the counter-term

```mathematica
amp1[0] = FCFAConvert[CreateFeynAmp[diag1[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {Zm -> SMP["Z_m"], Zphi -> SMP["Z_phi"], 
         GaugeXi[S[1]] -> 1, Mphi -> m}]
```

![1u8t2twsaydkj](img/1u8t2twsaydkj.svg)

Quartic vertex including the counter-term

```mathematica
amp2[0] = FCFAConvert[CreateFeynAmp[diag2[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p1}, OutgoingMomenta -> {p2, p3}, 
     LorentzIndexNames -> {mu}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {Zg -> SMP["Z_g"], Zphi -> SMP["Z_phi"], 
         GaugeXi[S[1]] -> 1, Mphi -> m}]
```

![1f00dlj1pfs95](img/1f00dlj1pfs95.svg)

## Calculate the amplitudes

### Self-energy

```mathematica
amp1[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 /. {SMP["Z_phi"] -> 1 + alpha*SMP["d_phi"], 
                  SMP["Z_m"] -> 1 + alpha*SMP["d_m"]} & )[amp1[0]]]]]
```

![0xxgvzjkoabqh](img/0xxgvzjkoabqh.svg)

Express the self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp1[2] = ToPaVe[amp1[1], l]
```

![1vmgidma16r76](img/1vmgidma16r76.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp1Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_m"], 
              SMP["d_phi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp1[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![06z2kzk565mtw](img/06z2kzk565mtw.svg)

Equating the result to zero and solving for d_phi and d_m we obtain  the renormalization constants in  the minimal subtraction schemes.

```mathematica
sol[1] = Simplify[Flatten[Solve[SelectNotFree2[amp1Div[0], 
               p] == 0, SMP["d_phi"]]]]; 
sol[2] = Simplify[Flatten[Solve[SelectFree2[amp1Div[0], p] == 
               0 /. sol[1], SMP["d_m"]]]]; 
solMS1 = Join[sol[1], sol[2]] /. 
     {SMP["d_phi"] -> SMP["d_phi^MS"], SMP["d_m"] -> 
         SMP["d_m^MS"], SMP["Delta"] -> 1/Epsilon}
solMSbar1 = Join[sol[1], sol[2]] /. 
     {SMP["d_phi"] -> SMP["d_phi^MSbar"], 
       SMP["d_m"] -> SMP["d_m^MSbar"]}
```

![0zwwun8gtx8w8](img/0zwwun8gtx8w8.svg)

![10riwqb40izsd](img/10riwqb40izsd.svg)

### Cubic vertex

```mathematica
amp2[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 //. {SMP["Z_g"] -> 1 + alpha*SMP["d_g"], 
                  SMP["Z_phi"] -> 1 + alpha*SMP["d_phi"]} & )[
     amp2[0]]]]]
```

![09igcnuipp6yz](img/09igcnuipp6yz.svg)

Express the cubic vertex in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp2[2] = ToPaVe[amp2[1], l]
```

![16zx9909ufnsp](img/16zx9909ufnsp.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp2Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_g"], 
              SMP["d_phi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp2[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![1aux2a1x25mgu](img/1aux2a1x25mgu.svg)

```mathematica
sol[3] = Simplify[Flatten[Solve[amp2Div[0] == 0 /. sol[1], 
           SMP["d_g"]]]]; 
solMS2 = sol[3] /. {SMP["d_g"] -> SMP["d_g^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar2 = sol[3] /. {SMP["d_g"] -> SMP["d_g^MSbar"]}
```

![0jrr2epyag3s8](img/0jrr2epyag3s8.svg)

![05erioxd1eo3d](img/05erioxd1eo3d.svg)

## Check the final results

```mathematica
knownResult = {SMP["d_phi^MS"] -> 0, SMP["d_m^MS"] -> 
         (g^2*(1/Epsilon))/(32*Pi^2*m^2), SMP["d_g^MS"] -> 0, 
       SMP["d_phi^MSbar"] -> 0, SMP["d_m^MSbar"] -> 
         (g^2*SMP["Delta"])/(32*Pi^2*m^2), SMP["d_g^MSbar"] -> 0}; 
FCCompareResults[Join[solMS1, solMS2, solMSbar1, solMSbar2], 
     knownResult, Text -> {"\tCompare to Cheng and Li, Gauge \
    theory of elementary particle physics, Problems and Solutions, \
    Eq. 2.120:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0sw3vrc99ulep](img/0sw3vrc99ulep.svg)

![03bnu8iyh2epa](img/03bnu8iyh2epa.svg)