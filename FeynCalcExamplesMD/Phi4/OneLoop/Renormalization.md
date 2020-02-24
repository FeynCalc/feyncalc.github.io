---
title: 1-loop phi^4 renormalization in the minimal subtraction schemes
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom phi^4 model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/Phi4/GenerateModelPhi4.m before running it for the first time.

```mathematica
description = "Renormalization, phi^4, MS and MSbar, 1-loop"; 
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
       Model -> FileNameJoin[{"Phi4", "Phi4"}], 
       GenericModel -> FileNameJoin[{"Phi4", "Phi4"}]}; 
top[i_, j_] := CreateTopologies[1, i -> j]; 
topCT[i_, j_] := CreateCTTopologies[1, i -> j]; 
topVertex[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> {WFCorrections}]; 
topVertexCT[i_, j_] := CreateCTTopologies[1, i -> j, 
       ExcludeTopologies -> {WFCorrectionCTs}]; 
{diagPhi4SE, diagPhi4SECT} = 
     (InsertFields[#1, {S[1]} -> {S[1]}, Sequence @@ 
              params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagVertex, diagVertexCT} = 
     (InsertFields[#1, {S[1], S[1]} -> {S[1], S[1]}, 
            Sequence @@ params] & ) /@ {topVertex[2, 2], 
         topVertexCT[2, 2]}; 
```

```mathematica
diag1[0] = diagPhi4SE[[0]][Sequence @@ diagPhi4SE, 
       Sequence @@ diagPhi4SECT]; 
diag2[0] = diagVertex[[0]][Sequence @@ diagVertex, 
       Sequence @@ diagVertexCT]; 
```

```mathematica
Paint[diag1[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![1gzkanb8o9vei](img/1gzkanb8o9vei.svg)

```mathematica
Paint[diag2[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![01ilidscsgyw2](img/01ilidscsgyw2.svg)

![01ffinsreq4eo](img/01ffinsreq4eo.svg)

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

![137kenedu0hw9](img/137kenedu0hw9.svg)

Quartic vertex including the counter-term

```mathematica
amp2[0] = FCFAConvert[CreateFeynAmp[diag2[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {p3, p4}, 
     LorentzIndexNames -> {mu}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {Zg -> SMP["Z_g"], Zphi -> SMP["Z_phi"], 
         GaugeXi[S[1]] -> 1, Mphi -> m}]
```

![1rnum8ru6stup](img/1rnum8ru6stup.svg)

## Calculate the amplitudes

### Self-energy

```mathematica
amp1[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 /. {SMP["Z_phi"] -> 1 + alpha*SMP["d_phi"], 
                  SMP["Z_m"] -> 1 + alpha*SMP["d_m"]} & )[amp1[0]]]]]
```

![0ajhf97r9fzbh](img/0ajhf97r9fzbh.svg)

Express the self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp1[2] = ToPaVe[amp1[1], l]
```

![0eqq2dsvd8gf1](img/0eqq2dsvd8gf1.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp1Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_m"], 
              SMP["d_phi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp1[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![08qhr3s98rp5e](img/08qhr3s98rp5e.svg)

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

![17ecpy4rv03o2](img/17ecpy4rv03o2.svg)

![1oir2o54k0btz](img/1oir2o54k0btz.svg)

### Quartic vertex

```mathematica
amp2[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 //. {SMP["Z_g"] -> 1 + alpha*SMP["d_g"], 
                  SMP["Z_phi"] -> 1 + alpha*SMP["d_phi"]} & )[
     amp2[0]]]]]
```

![17c8bdp46gfjq](img/17c8bdp46gfjq.svg)

Express the quartic vertex in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp2[2] = ToPaVe[amp2[1], l]
```

![02g3iqedm1cb5](img/02g3iqedm1cb5.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp2Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_g"], 
              SMP["d_phi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp2[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![0s4dvghtq1xui](img/0s4dvghtq1xui.svg)

```mathematica
sol[3] = Simplify[Flatten[Solve[amp2Div[0] == 0 /. sol[1], 
           SMP["d_g"]]]]; 
solMS2 = sol[3] /. {SMP["d_g"] -> SMP["d_g^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar2 = sol[3] /. {SMP["d_g"] -> SMP["d_g^MSbar"]}
```

![1j7ow5ei3qemp](img/1j7ow5ei3qemp.svg)

![0ghjzw4tex5iq](img/0ghjzw4tex5iq.svg)

## Check the final results

```mathematica
knownResult = {SMP["d_phi^MS"] -> 0, SMP["d_m^MS"] -> 
         (g*(1/Epsilon))/(32*Pi^2), SMP["d_g^MS"] -> 
         (3*g*(1/Epsilon))/(32*Pi^2), SMP["d_phi^MSbar"] -> 0, 
       SMP["d_m^MSbar"] -> (g*SMP["Delta"])/(32*Pi^2), 
       SMP["d_g^MSbar"] -> (3*g*SMP["Delta"])/(32*Pi^2)}; 
FCCompareResults[Join[solMS1, solMS2, solMSbar1, solMSbar2], 
     knownResult, Text -> {"\tCompare to Bailin and Love, \
    Introduction to Gauge Field Theory, Eqs. 7.73-7.74 and Eqs. \
    7.76-7.77:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1svptvx4xab2x](img/1svptvx4xab2x.svg)

![0ez0hkopy5yua](img/0ez0hkopy5yua.svg)