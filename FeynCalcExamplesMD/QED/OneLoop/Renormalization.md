---
title: 1-loop QED renormalization in the minimal subtraction schemes
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom QED model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/QED/GenerateModelQED.m before running it for the first time.

```mathematica
description = "Renormalization, QED, MS and MSbar, 1-loop"; 
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

Nicer typesetting

```mathematica
MakeBoxes[mu, TraditionalForm] := "\[Mu]"; 
MakeBoxes[nu, TraditionalForm] := "\[Nu]"; 
```

```mathematica
params = {InsertionLevel -> {Particles}, 
       Model -> FileNameJoin[{"QED", "QED"}], 
       GenericModel -> FileNameJoin[{"QED", "QED"}], 
       ExcludeParticles -> {F[2, {2 | 3}]}}; 
top[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs}]; 
topCT[i_, j_] := CreateCTTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs}]; 
{diagElectronSE, diagElectronSECT} = 
     (InsertFields[#1, {F[2, {1}]} -> {F[2, {1}]}, 
            Sequence @@ params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagPhotonSE, diagPhotonSECT} = 
     (InsertFields[#1, {V[1]} -> {V[1]}, Sequence @@ 
              params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagVertex, diagVertexCT} = 
     (InsertFields[#1, {F[2, {1}], V[1]} -> {F[2, {1}]}, 
            Sequence @@ params] & ) /@ {top[2, 1], topCT[2, 1]}; 
```

```mathematica
diag1[0] = diagElectronSE[[0]][diagElectronSE[[1]], 
       diagElectronSECT[[1]]]; 
diag2[0] = diagPhotonSE[[0]][diagPhotonSE[[1]], 
       diagPhotonSECT[[1]]]; 
diag3[0] = diagVertex[[0]][diagVertex[[1]], diagVertexCT[[1]]]; 
```

```mathematica
Paint[diag1[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![1juen4fabkjsg](img/1juen4fabkjsg.svg)

```mathematica
Paint[diag2[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![0msf65weytixb](img/0msf65weytixb.svg)

```mathematica
Paint[diag3[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![1ubkz9xbago5e](img/1ubkz9xbago5e.svg)

## Obtain the amplitudes

The 1/(2Pi)^D prefactor is implicit. We need to replace e with -e to be compatible
with the convention D^mu = d^mu + ie A^mu

Electron self-energy including the counter-term

```mathematica
amp1[0] = FCFAConvert[CreateFeynAmp[diag1[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {Zm -> SMP["Z_m"], Zpsi -> SMP["Z_psi"], 
         SMP["e"] -> Sqrt[4*Pi*SMP["alpha_fs"]], 
         GaugeXi[V[1]] -> GaugeXi}, Contract -> True]
```

![10k2y99rpjxsn](img/10k2y99rpjxsn.svg)

Photon self-energy including the counter-term

```mathematica
amp2[0] = FCTraceFactor[FCFAConvert[CreateFeynAmp[diag2[0], 
         Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
       IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
       LorentzIndexNames -> {mu, nu}, LoopMomenta -> {l}, 
       UndoChiralSplittings -> True, ChangeDimension -> D, 
       List -> False, SMP -> True, FinalSubstitutions -> 
         {ZA -> SMP["Z_A"], Zxi -> SMP["Z_xi"], 
           SMP["e"] -> Sqrt[4*Pi*SMP["alpha_fs"]], 
           GaugeXi[V[1]] -> GaugeXi}, Contract -> True]]
```

![1ea5z3lmx1j02](img/1ea5z3lmx1j02.svg)

Electron-photon vertex including the counter-term

```mathematica
amp3[0] = FCFAConvert[CreateFeynAmp[diag3[0], 
         Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
       IncomingMomenta -> {p1, k}, OutgoingMomenta -> {p2}, 
       LorentzIndexNames -> {mu}, LoopMomenta -> {l}, 
       UndoChiralSplittings -> True, ChangeDimension -> D, 
       List -> False, SMP -> True, FinalSubstitutions -> 
         {ZA -> SMP["Z_A"], Ze -> SMP["Z_e"], Zpsi -> SMP["Z_psi"], 
           SMP["e"]^3 -> 4*Pi*SMP["alpha_fs"]*SMP["e"], 
           GaugeXi[V[1]] -> GaugeXi}, Contract -> True] /. 
     SMP["e"] -> -SMP["e"]
```

![1b9x7hn8h5e3w](img/1b9x7hn8h5e3w.svg)

## Calculate the amplitudes

### Electron self-energy

```mathematica
amp1[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 /. {SMP["Z_psi"] -> 1 + alpha*SMP["d_psi"], 
                  SMP["Z_m"] -> 1 + alpha*SMP["d_m"]} & )[amp1[0]]]]]
```

![08gpktmyysi6b](img/08gpktmyysi6b.svg)

Tensor reduction allows us to express the electron self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp1[2] = TID[amp1[1], l, ToPaVe -> True]
```

![0xgecgl01cpg6](img/0xgecgl01cpg6.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp1Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_m"], 
              SMP["d_psi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp1[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![0q8cqkutbhotw](img/0q8cqkutbhotw.svg)

Equating the result to zero and solving for d_psi and d_m we obtain  the renormalization constants in 
the minimal subtraction schemes.

```mathematica
sol[1] = Simplify[Flatten[Solve[SelectNotFree2[amp1Div[0], 
               DiracGamma] == 0, SMP["d_psi"]]]]; 
sol[2] = Simplify[Flatten[
         Solve[SelectFree2[amp1Div[0], DiracGamma] == 0 /. sol[1], 
           SMP["d_m"]]]]; 
solMS1 = Join[sol[1], sol[2]] /. 
     {SMP["d_psi"] -> SMP["d_psi^MS"], SMP["d_m"] -> 
         SMP["d_m^MS"], SMP["Delta"] -> 1/Epsilon}
solMSbar1 = Join[sol[1], sol[2]] /. 
     {SMP["d_psi"] -> SMP["d_psi^MSbar"], 
       SMP["d_m"] -> SMP["d_m^MSbar"]}
```

![02xv8klym31eo](img/02xv8klym31eo.svg)

![0glkim02hhmop](img/0glkim02hhmop.svg)

### Photon self-energy

```mathematica
amp2[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 //. {SMP["Z_xi"] -> SMP["Z_A"], SMP["Z_A"] -> 
                    1 + alpha*SMP["d_A"]} & )[amp2[0]]]]]
```

![0qvhy8ehv7wpa](img/0qvhy8ehv7wpa.svg)

Tensor reduction allows us to express the electron self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
amp2[2] = TID[amp2[1], l, ToPaVe -> True]
```

![0lvqnqpmn3z6h](img/0lvqnqpmn3z6h.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp2Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_A"]}] & )[
       FCHideEpsilon[Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               PaVeUVPart[amp2[2], Prefactor -> 1/(2*Pi)^D]]]]]]]
```

![1qsqglhwc2fuk](img/1qsqglhwc2fuk.svg)

Equating this to zero and solving for d_A we obtain the wave-function renormalization constant for the photon in the minimal subtraction schemes.

```mathematica
sol[3] = Flatten[Solve[amp2Div[0] == 0, SMP["d_A"]]]; 
solMS2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MSbar"]}
```

![0lybvk09ukn4s](img/0lybvk09ukn4s.svg)

![0iyp6ka7fry53](img/0iyp6ka7fry53.svg)

### Electron-photon vertex

```mathematica
amp3[1] = (#1 /. alpha -> 1 & )[
     Normal[(Series[#1, {alpha, 0, 1}] & )[
         (#1 //. {SMP["Z_psi"] -> 1 + alpha*SMP["d_psi"], 
                  SMP["Z_A"] -> 1 + alpha*SMP["d_A"], SMP["Z_e"] -> 
                    1 + alpha*SMP["d_e"]} & )[amp3[0]]]]]
```

![0ddyphusm37rm](img/0ddyphusm37rm.svg)

The result of the tensor reduction is quite large, since we keep the full gauge dependence and do not specify the kinematics

```mathematica
amp3[2] = TID[amp3[1], l, ToPaVe -> True, UsePaVeBasis -> True]
```

![1x8u7i7dundek](img/1x8u7i7dundek.svg)

Discard all the finite pieces of the 1-loop amplitude

```mathematica
amp3Div[0] = Simplify[
     (SelectNotFree2[#1, {SMP["Delta"], SMP["d_A"], SMP["d_e"], 
              SMP["d_psi"]}] & )[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
               DiracSimplify[PaVeUVPart[amp3[2], Prefactor -> 
                     1/(2*Pi)^D]]]]]]]]
```

![1g39jxov25kyo](img/1g39jxov25kyo.svg)

```mathematica
ward[0] = Simplify[amp3Div[0]/(-FCI[I*SMP["e"]*GAD[mu]]) == 0]
```

![03s93njbiddvh](img/03s93njbiddvh.svg)

```mathematica
wardMS[0] = Simplify[ward[0] /. Epsilon -> 1/SMP["Delta"] /. 
         {SMP["d_psi"] -> SMP["d_psi^MSbar"]} /. solMSbar1]
wardMSbar[0] = Simplify[
     ward[0] /. {SMP["d_psi"] -> SMP["d_psi^MSbar"]} /. solMSbar1]
```

![1iiupc2y88f2k](img/1iiupc2y88f2k.svg)

![1kfmp5pbe1eqh](img/1kfmp5pbe1eqh.svg)

```mathematica
knownResults = {SMP["d_A"] + 2*SMP["d_e"] == 0, 
       SMP["d_A"] + 2*SMP["d_e"] == 0}; 
FCCompareResults[{wardMS[0], wardMSbar[0]}, knownResults, 
     Text -> {"\tVerify Ward's identity:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
```

![053m9cxkra8ie](img/053m9cxkra8ie.svg)

## Check the final results

```mathematica
knownResult = {SMP["d_psi^MS"] -> -(GaugeXi*SMP["alpha_fs"])/
           (4*Epsilon*Pi), SMP["d_m^MS"] -> (-3*SMP["alpha_fs"])/
           (4*Epsilon*Pi), SMP["d_A^MS"] -> -SMP["alpha_fs"]/
           (3*Epsilon*Pi), SMP["d_psi^MSbar"] -> 
         -(GaugeXi*SMP["alpha_fs"]*SMP["Delta"])/(4*Pi), 
       SMP["d_m^MSbar"] -> (-3*SMP["alpha_fs"]*SMP["Delta"])/
           (4*Pi), SMP["d_A^MSbar"] -> 
         -(SMP["alpha_fs"]*SMP["Delta"])/(3*Pi)}; 
FCCompareResults[Join[solMS1, solMS2, solMSbar1, solMSbar2], 
     knownResult, Text -> {"\tCheck the final result:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0mv6mshl1ar9c](img/0mv6mshl1ar9c.svg)

![0i0ff662q8xe7](img/0i0ff662q8xe7.svg)