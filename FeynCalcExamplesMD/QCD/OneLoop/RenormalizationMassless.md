---
title: 1-loop massless QCD renormalization in the minimal subtraction schemes
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom QCD model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/QCD/GenerateModelQCD.m before running it for the first time.

```mathematica
description = 
     "Renormalization, massless QCD, MS and MSbar, 1-loop"; 
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
MakeBoxes[rho, TraditionalForm] := "\[Rho]"; 
MakeBoxes[si, TraditionalForm] := "\[Sigma]"; 
```

```mathematica
params = {InsertionLevel -> {Particles}, 
       Model -> FileNameJoin[{"QCD", "QCD"}], 
       GenericModel -> FileNameJoin[{"QCD", "QCD"}], 
       ExcludeParticles -> {F[3 | 4, {2 | 3}], F[4, {1}]}}; 
top[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs}]; 
topTriangle[i_, j_] := CreateTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs, SelfEnergies}]; 
topCT[i_, j_] := CreateCTTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs}]; 
topTriangleCT[i_, j_] := CreateCTTopologies[1, i -> j, 
       ExcludeTopologies -> {Tadpoles, WFCorrections, 
           WFCorrectionCTs, SelfEnergyCTs}]; 
{diagQuarkSE, diagQuarkSECT} = 
     (InsertFields[#1, {F[3, {1}]} -> {F[3, {1}]}, 
            Sequence @@ params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagGluonSE, diagGluonSECT} = 
     (InsertFields[#1, {V[5]} -> {V[5]}, Sequence @@ 
              params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagGhostSE, diagGhostSECT} = 
     (InsertFields[#1, {U[5]} -> {U[5]}, Sequence @@ 
              params] & ) /@ {top[1, 1], topCT[1, 1]}; 
{diagQuarkGluonVertex, diagQuarkGluonVertexCT} = 
     (InsertFields[#1, {F[3, {1}], V[5]} -> {F[3, {1}]}, 
            Sequence @@ params] & ) /@ {topTriangle[2, 1], 
         topTriangleCT[2, 1]}; 
```

```mathematica
diag1[0] = diagQuarkSE[[0]][Sequence @@ diagQuarkSE, 
       Sequence @@ diagQuarkSECT]; 
diag2[0] = diagGluonSE[[0]][Sequence @@ diagGluonSE, 
       Sequence @@ diagGluonSECT]; 
diag3[0] = diagGhostSE[[0]][Sequence @@ diagGhostSE, 
       Sequence @@ diagGhostSECT]; 
diag4[0] = diagQuarkGluonVertex[[0]][
       Sequence @@ diagQuarkGluonVertex, 
       Sequence @@ diagQuarkGluonVertexCT]; 
```

```mathematica
Paint[diag1[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![02anaixr2dhbv](img/02anaixr2dhbv.svg)

```mathematica
Paint[diag2[0], ColumnsXRows -> {4, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 128}]; 
```

![1xoxiz04sycl9](img/1xoxiz04sycl9.svg)

![1upornmo194rx](img/1upornmo194rx.svg)

```mathematica
Paint[diag3[0], ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![0kinr7l36muip](img/0kinr7l36muip.svg)

```mathematica
Paint[diag4[0], ColumnsXRows -> {3, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![09pnf715m5lm7](img/09pnf715m5lm7.svg)

## Obtain the amplitudes

The 1/(2Pi)^D prefactor is implicit.

Quark self-energy including the counter-term

```mathematica
ampQuarkSE[0] = FCFAConvert[CreateFeynAmp[diag1[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu, nu}, DropSumOver -> True, 
     LoopMomenta -> {l}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> False, SMP -> True, 
     FinalSubstitutions -> {Zm -> SMP["Z_m"], 
         Zpsi -> SMP["Z_psi"], SMP["m_u"] -> 0}]
```

![1xdbt78qm43sh](img/1xdbt78qm43sh.svg)

Gluon self-energy including the counter-term

```mathematica
ampGluonSE[0] = FCFAConvert[CreateFeynAmp[diag2[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu, nu, rho, si}, DropSumOver -> True, 
     LoopMomenta -> {l}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> True, SMP -> True, 
     FinalSubstitutions -> {ZA -> SMP["Z_A"], Zxi -> SMP["Z_xi"], 
         SMP["m_u"] -> 0}]
```

![04xximvzm4p7r](img/04xximvzm4p7r.svg)

Ghost self-energy including the counter-term

```mathematica
ampGhostSE[0] = FCFAConvert[CreateFeynAmp[diag3[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu, nu}, DropSumOver -> True, 
     LoopMomenta -> {l}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> False, SMP -> True, 
     FinalSubstitutions -> {Zu -> SMP["Z_u"]}]
```

![0z103j1i949y7](img/0z103j1i949y7.svg)

Quark-gluon vertex including the counter-term

```mathematica
ampQGlVertex[0] = FCFAConvert[CreateFeynAmp[diag4[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p1, k}, OutgoingMomenta -> {p2}, 
     LorentzIndexNames -> {mu, nu, rho}, DropSumOver -> True, 
     LoopMomenta -> {l}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> False, SMP -> True, 
     FinalSubstitutions -> {ZA -> SMP["Z_A"], Zg -> SMP["Z_g"], 
         Zpsi -> SMP["Z_psi"], SMP["m_u"] -> 0}]
```

![0cqpfruwix5oi](img/0cqpfruwix5oi.svg)

## Calculate the amplitudes

### Quark self-energy

Tensor reduction allows us to express the quark self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
ampQuarkSE[1] = (TID[#1, l, UsePaVeBasis -> True, 
            ToPaVe -> True] & )[DiracSimplify[
         SUNSimplify[ampQuarkSE[0]]]]; 
```

Discard all the finite pieces of the 1-loop amplitude.

```mathematica
ampQuarkSEDiv[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & )[
       ampQuarkSE[1]]; 
```

```mathematica
ampQuarkSEDiv[1] = Simplify[FCHideEpsilon[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           FCReplaceD[ampQuarkSEDiv[0], D -> 4 - 2*Epsilon]]]]]
```

![0axvwbjppiq72](img/0axvwbjppiq72.svg)

```mathematica
ampQuarkSEDiv[2] = (SelectNotFree2[#1, SMP["Delta"], 
          SMP["d_m"], SMP["d_psi"]] & )[(#1 /. alpha -> 1 & )[
       Normal[(Series[#1, {alpha, 0, 1}] & )[
           (#1 //. {SMP["Z_m"] -> 1 + alpha*SMP["d_m"], 
                    SMP["Z_psi"] -> 1 + alpha*SMP["d_psi"]} & )[
             ampQuarkSEDiv[1]]]]]]
```

![0ku3x29s7a8xi](img/0ku3x29s7a8xi.svg)

```mathematica
ampQuarkSEDiv[3] = 
   (Collect2[#1, DiracGamma, Factoring -> Simplify] & )[
     SUNSimplify[ampQuarkSEDiv[2]]]
```

![0mlxx2rvu347z](img/0mlxx2rvu347z.svg)

```mathematica
sol[1] = (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       (#1 /. (a_ -> b_) :> a -> SUNSimplify[b] & )[
         Flatten[Solve[ampQuarkSEDiv[3] == 0, SMP["d_psi"]]]]]; 
solMS1 = sol[1] /. {SMP["d_psi"] -> SMP["d_psi^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar1 = sol[1] /. {SMP["d_psi"] -> SMP["d_psi^MSbar"]}
```

![11as8s6vaalub](img/11as8s6vaalub.svg)

![01mtod2cakon8](img/01mtod2cakon8.svg)

### Gluon self-energy

Tensor reduction allows us to express the gluon self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
ampGluonSE[1] = DiracSimplify[SUNSimplify[ampGluonSE[0][[1]] + 
           Nf*ampGluonSE[0][[2]] + Total[ampGluonSE[0][[
                3 ;; All]]]]]; 
```

```mathematica
ampGluonSE[2] = TID[ampGluonSE[1], l, UsePaVeBasis -> True, 
       ToPaVe -> True]; 
```

Discard all the finite pieces of the 1-loop amplitude

```mathematica
ampGluonSEDiv[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & )[
     ampGluonSE[2]]
```

![0rwkbi86mv962](img/0rwkbi86mv962.svg)

```mathematica
ampGluonSEDiv[1] = SUNSimplify[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             FCReplaceD[ampGluonSEDiv[0], D -> 4 - 2*Epsilon]]]]]; 
```

```mathematica
ampGluonSEDiv[2] = (SelectNotFree2[#1, SMP["Delta"], 
          SMP["d_A"], SMP["d_xi"]] & )[(#1 /. alpha -> 1 & )[
       Normal[(Series[#1, {alpha, 0, 1}] & )[
           (#1 //. {SMP["Z_A"] -> 1 + alpha*SMP["d_A"], 
                    SMP["Z_xi"] -> 1 + alpha*SMP["d_A"]} & )[
             ampGluonSEDiv[1]]]]]]
```

![0xwjp7e9qb9pd](img/0xwjp7e9qb9pd.svg)

```mathematica
sol[3] = Simplify[
     (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       Flatten[Solve[ampGluonSEDiv[2] == 0, SMP["d_A"]]]]]
solMS2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MSbar"]}
```

![1y1ems6rhrmon](img/1y1ems6rhrmon.svg)

![1riitq7qlga1q](img/1riitq7qlga1q.svg)

![0ixr19wc6zhvt](img/0ixr19wc6zhvt.svg)

### Ghost self-energy

Tensor reduction allows us to express the ghost self-energy in tems of the Passarino-Veltman coefficient functions.

```mathematica
ampGhostSE[1] = DiracSimplify[SUNSimplify[ampGhostSE[0]]]; 
```

```mathematica
ampGhostSE[2] = TID[ampGhostSE[1], l, UsePaVeBasis -> True, 
       ToPaVe -> True]; 
```

Discard all the finite pieces of the 1-loop amplitude

```mathematica
ampGhostSEDiv[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & )[
     ampGhostSE[2]]
```

![015f87h6c4cd4](img/015f87h6c4cd4.svg)

```mathematica
ampGhostSEDiv[1] = SUNSimplify[FCHideEpsilon[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           FCReplaceD[ampGhostSEDiv[0], D -> 4 - 2*Epsilon]]]]]
```

![1be44ws01vi9x](img/1be44ws01vi9x.svg)

```mathematica
ampGhostSEDiv[2] = Simplify[
     (SelectNotFree2[#1, SMP["Delta"], SMP["d_u"]] & )[
       (#1 /. alpha -> 1 & )[Normal[(Series[#1, {alpha, 0, 1}] & )[
             (#1 //. {SMP["Z_u"] -> 1 + alpha*SMP["d_u"]} & )[
               ampGhostSEDiv[1]]]]]]]
```

![064cjlkvvw1bm](img/064cjlkvvw1bm.svg)

```mathematica
sol[4] = Simplify[
     (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       Flatten[Solve[ampGhostSEDiv[2] == 0, SMP["d_u"]]]]]
solMS3 = sol[4] /. {SMP["d_u"] -> SMP["d_u^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar3 = sol[4] /. {SMP["d_u"] -> SMP["d_u^MSbar"]}
```

![06hqis3cojmht](img/06hqis3cojmht.svg)

![0vepwbvjnxpvb](img/0vepwbvjnxpvb.svg)

![0hhxpl6gcx264](img/0hhxpl6gcx264.svg)

### Quark-gluon vertex

Tensor reduction allows us to express the quark-gluon vertex in tems of the Passarino-Veltman coefficient functions.

```mathematica
ampQGlVertex[1] = DiracSimplify[SUNSimplify[ampQGlVertex[0]]]; 
```

```mathematica
ampQGlVertex[2] = TID[ampQGlVertex[1], l, UsePaVeBasis -> True, 
       ToPaVe -> True]; 
```

Discard all the finite pieces of the 1-loop amplitude

```mathematica
ampQGlVertexDiv[0] = 
   (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & )[ampQGlVertex[2]]
```

![1x2hvdtg7eb2r](img/1x2hvdtg7eb2r.svg)

```mathematica
ampQGlVertexDiv[1] = SUNSimplify[FCHideEpsilon[
         Normal[(Series[#1, {Epsilon, 0, 0}] & )[
             FCReplaceD[ampQGlVertexDiv[0], D -> 4 - 2*Epsilon]]]]]; 
```

```mathematica
ampQGlVertexDiv[2] = Simplify[
     (SelectNotFree2[#1, SMP["Delta"], SMP["d_g"], SMP["d_A"], 
            SMP["d_psi"]] & )[(#1 /. alpha -> 1 & )[
         Normal[(Series[#1, {alpha, 0, 1}] & )[
             (#1 //. {SMP["Z_g"] -> 1 + alpha*SMP["d_g"], 
                      SMP["Z_A"] -> 1 + alpha*SMP["d_A"], 
           SMP["Z_psi"] -> 
                        1 + alpha*SMP["d_psi"]} & )[
       ampQGlVertexDiv[1]]]]]]]
```

![07py4rlh4pj7i](img/07py4rlh4pj7i.svg)

```mathematica
ampQGlVertexDiv[3] = (Collect2[#1, Epsilon, SUNIndex] & )[
     (#1 /. SUNTrace[x__] :> SUNTrace[x, Explicit -> True] & )[
       (SUNSimplify[#1, Explicit -> True] & )[ampQGlVertexDiv[2]]]]
```

![1hk31drm2g7ki](img/1hk31drm2g7ki.svg)

```mathematica
ampQGlVertexDiv[4] = (Collect2[#1, SMP] & )[
     SUNSimplify[(#1 /. suntf[xx_, _SUNFIndex, _SUNFIndex] :> 
                SUNT @@ xx & )[(#1 /. SUNTF -> suntf & )[
           (#1 /. suntf[xx_, _SUNFIndex, _SUNFIndex] :> 
                    SUNT @@ xx & )[ampQGlVertexDiv[3]]]]]]
```

![1uwpnwjgjyjef](img/1uwpnwjgjyjef.svg)

```mathematica
ampQGlVertexDiv[5] = SUNSimplify[(Collect2[#1, Epsilon] & )[
       (#1 /. SMP["g_s"]^3 -> 4*Pi*SMP["alpha_s"]*SMP["g_s"] & )[
         ampQGlVertexDiv[4] /. {SMP["d_A"] -> SMP["d_A^MS"], 
                 SMP["d_psi"] -> SMP["d_psi^MS"], SMP["d_g"] -> 
                   SMP["d_g"], SMP["Delta"] -> 1/Epsilon} /. 
      solMS1 /. 
           solMS2]]]
```

![04r62uk8fs9pt](img/04r62uk8fs9pt.svg)

```mathematica
sol[5] = Simplify[Flatten[Solve[ampQGlVertexDiv[5] == 0, 
         SMP["d_g"]]]]
solMS4 = sol[5] /. {SMP["d_g"] -> SMP["d_g^MS"]}
solMSbar4 = sol[5] /. {SMP["d_g"] -> SMP["d_g^MSbar"], 
       1/Epsilon -> SMP["Delta"]}
```

![0gmitdm7cjw01](img/0gmitdm7cjw01.svg)

![0op2fjwest78h](img/0op2fjwest78h.svg)

![1334jdt3nbz6q](img/1334jdt3nbz6q.svg)

## Check the final results

```mathematica
knownResult = Factor2[{SMP["d_psi^MS"] -> 
           (-SMP["alpha_s"]/(4*Pi))*(1/Epsilon)*CF*GaugeXi["G"], 
         SMP["d_psi^MSbar"] -> (-SMP["alpha_s"]/(4*Pi))*
             SMP["Delta"]*CF*GaugeXi["G"], SMP["d_A^MS"] -> 
           (SMP["alpha_s"]/(4*Pi))*(1/Epsilon)*
             ((1/2)*CA*(13/3 - GaugeXi["G"]) - (2/3)*Nf), 
         SMP["d_A^MSbar"] -> (SMP["alpha_s"]/(4*Pi))*SMP["Delta"]*
             ((1/2)*CA*(13/3 - GaugeXi["G"]) - (2/3)*Nf), 
         SMP["d_u^MS"] -> (SMP["alpha_s"]/(4*Pi))*(CA/Epsilon)*
             ((3 - GaugeXi["G"])/4), SMP["d_u^MSbar"] -> 
           (SMP["alpha_s"]/(4*Pi))*CA*SMP["Delta"]*
             ((3 - GaugeXi["G"])/4), SMP["d_g^MS"] -> 
           (-11*CA*SMP["alpha_s"])/(24*Epsilon*Pi) + 
             (Nf*SMP["alpha_s"])/(12*Epsilon*Pi), 
         SMP["d_g^MSbar"] -> ((-11*CA*SMP["alpha_s"])/(24*Pi))*
               SMP["Delta"] + ((Nf*SMP["alpha_s"])/(12*Pi))*
               SMP["Delta"]}]; 
FCCompareResults[Factor2[Join[solMS1, solMSbar1, solMS2, 
       solMSbar2, solMS3, solMSbar3, solMS4, solMSbar4]], 
   knownResult, Text -> 
     {"\tCompare to Muta, Foundations of QCD, Eqs \
   2.5.131-2.5.147:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1osomhsxziqf6](img/1osomhsxziqf6.svg)

![05il09pa8ck2w](img/05il09pa8ck2w.svg)

![1fmut8ap4815k](img/1fmut8ap4815k.svg)