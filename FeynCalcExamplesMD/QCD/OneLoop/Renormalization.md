---
title: 1-loop QCD renormalization in the minimal subtraction schemes
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom QCD model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/QCD/GenerateModelQCD.m before running it for the first time.

```mathematica
description = "Renormalization, QCD, MS and MSbar, 1-loop"; 
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
         Zpsi -> SMP["Z_psi"], SMP["m_u"] -> SMP["m_q"]}]
```

![07sxtdygdb95j](img/07sxtdygdb95j.svg)

Gluon self-energy including the counter-term

```mathematica
ampGluonSE[0] = FCFAConvert[CreateFeynAmp[diag2[0], 
       Truncated -> True, GaugeRules -> {}, PreFactor -> 1], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {p}, 
     LorentzIndexNames -> {mu, nu, rho, si}, DropSumOver -> True, 
     LoopMomenta -> {l}, UndoChiralSplittings -> True, 
     ChangeDimension -> D, List -> True, SMP -> True, 
     FinalSubstitutions -> {ZA -> SMP["Z_A"], Zxi -> SMP["Z_xi"], 
         SMP["m_u"] -> SMP["m_q"]}]
```

![1agkdj2ar5tzu](img/1agkdj2ar5tzu.svg)

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
         Zpsi -> SMP["Z_psi"], SMP["m_u"] -> SMP["m_q"]}]
```

![1jecmxjjdlrjp](img/1jecmxjjdlrjp.svg)

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

![16bp8dldsb761](img/16bp8dldsb761.svg)

```mathematica
ampQuarkSEDiv[2] = (SelectNotFree2[#1, SMP["Delta"], 
          SMP["d_m"], SMP["d_psi"]] & )[(#1 /. alpha -> 1 & )[
       Normal[(Series[#1, {alpha, 0, 1}] & )[
           (#1 //. {SMP["Z_m"] -> 1 + alpha*SMP["d_m"], 
                    SMP["Z_psi"] -> 1 + alpha*SMP["d_psi"]} & )[
             ampQuarkSEDiv[1]]]]]]
```

![1v1nha7cztnyh](img/1v1nha7cztnyh.svg)

```mathematica
ampQuarkSEDiv[3] = 
   (Collect2[#1, DiracGamma, Factoring -> Simplify] & )[
     SUNSimplify[ampQuarkSEDiv[2]]]
```

![0br64x3fksrrh](img/0br64x3fksrrh.svg)

```mathematica
sol[1] = (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       (#1 /. (a_ -> b_) :> a -> SUNSimplify[b] & )[
         Flatten[Solve[SelectNotFree2[ampQuarkSEDiv[3], 
                 DiracGamma] == 0, SMP["d_psi"]]]]]; 
sol[2] = (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       (#1 /. (a_ -> b_) :> a -> SUNSimplify[b] & )[
         Flatten[Solve[SelectFree2[ampQuarkSEDiv[3], DiracGamma] == 
                 0 /. sol[1], SMP["d_m"]]]]]; 
solMS1 = Join[sol[1], sol[2]] /. 
     {SMP["d_psi"] -> SMP["d_psi^MS"], SMP["d_m"] -> 
         SMP["d_m^MS"], SMP["Delta"] -> 1/Epsilon}
solMSbar1 = Join[sol[1], sol[2]] /. 
     {SMP["d_psi"] -> SMP["d_psi^MSbar"], 
       SMP["d_m"] -> SMP["d_m^MSbar"]}
```

![07ckcafx5q21s](img/07ckcafx5q21s.svg)

![1cymlylg2b9uw](img/1cymlylg2b9uw.svg)

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

![1a7nj0tr538fx](img/1a7nj0tr538fx.svg)

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

![0fd0xpmeujy9p](img/0fd0xpmeujy9p.svg)

```mathematica
sol[3] = Simplify[
     (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       Flatten[Solve[ampGluonSEDiv[2] == 0, SMP["d_A"]]]]]
solMS2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar2 = sol[3] /. {SMP["d_A"] -> SMP["d_A^MSbar"]}
```

![0bufs8f19t8o8](img/0bufs8f19t8o8.svg)

![0t2p3ze5pcuzv](img/0t2p3ze5pcuzv.svg)

![0yufnw1dg93dr](img/0yufnw1dg93dr.svg)

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

![13doupd05a4w5](img/13doupd05a4w5.svg)

```mathematica
ampGhostSEDiv[1] = SUNSimplify[FCHideEpsilon[
       Normal[(Series[#1, {Epsilon, 0, 0}] & )[
           FCReplaceD[ampGhostSEDiv[0], D -> 4 - 2*Epsilon]]]]]
```

![1tj7e5295sofh](img/1tj7e5295sofh.svg)

```mathematica
ampGhostSEDiv[2] = Simplify[
     (SelectNotFree2[#1, SMP["Delta"], SMP["d_u"]] & )[
       (#1 /. alpha -> 1 & )[Normal[(Series[#1, {alpha, 0, 1}] & )[
             (#1 //. {SMP["Z_u"] -> 1 + alpha*SMP["d_u"]} & )[
               ampGhostSEDiv[1]]]]]]]
```

![1640jt1rhyv9f](img/1640jt1rhyv9f.svg)

```mathematica
sol[4] = Simplify[
     (#1 /. SMP["g_s"]^2 -> 4*Pi*SMP["alpha_s"] & )[
       Flatten[Solve[ampGhostSEDiv[2] == 0, SMP["d_u"]]]]]
solMS3 = sol[4] /. {SMP["d_u"] -> SMP["d_u^MS"], 
       SMP["Delta"] -> 1/Epsilon}
solMSbar3 = sol[4] /. {SMP["d_u"] -> SMP["d_u^MSbar"]}
```

![0gwnycpjs9a9z](img/0gwnycpjs9a9z.svg)

![12damdrv6m0fs](img/12damdrv6m0fs.svg)

![0worzq6dlrika](img/0worzq6dlrika.svg)

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

![1f0zw0db4jr25](img/1f0zw0db4jr25.svg)

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

![0c0njpjkw5o2c](img/0c0njpjkw5o2c.svg)

```mathematica
ampQGlVertexDiv[3] = (Collect2[#1, Epsilon, SUNIndex] & )[
     (#1 /. SUNTrace[x__] :> SUNTrace[x, Explicit -> True] & )[
       (SUNSimplify[#1, Explicit -> True] & )[ampQGlVertexDiv[2]]]]
```

![0jfhzvkr59g7a](img/0jfhzvkr59g7a.svg)

```mathematica
ampQGlVertexDiv[4] = (Collect2[#1, SMP] & )[
     SUNSimplify[(#1 /. suntf[xx_, _SUNFIndex, _SUNFIndex] :> 
                SUNT @@ xx & )[(#1 /. SUNTF -> suntf & )[
           (#1 /. suntf[xx_, _SUNFIndex, _SUNFIndex] :> 
                    SUNT @@ xx & )[ampQGlVertexDiv[3]]]]]]
```

![0nm52o3nyuo1e](img/0nm52o3nyuo1e.svg)

```mathematica
ampQGlVertexDiv[5] = SUNSimplify[(Collect2[#1, Epsilon] & )[
       (#1 /. SMP["g_s"]^3 -> 4*Pi*SMP["alpha_s"]*SMP["g_s"] & )[
         ampQGlVertexDiv[4] /. {SMP["d_A"] -> SMP["d_A^MS"], 
                 SMP["d_psi"] -> SMP["d_psi^MS"], SMP["d_g"] -> 
                   SMP["d_g"], SMP["Delta"] -> 1/Epsilon} /. 
      solMS1 /. 
           solMS2]]]
```

![162i1cr3ais07](img/162i1cr3ais07.svg)

```mathematica
sol[5] = Simplify[Flatten[Solve[ampQGlVertexDiv[5] == 0, 
         SMP["d_g"]]]]
solMS4 = sol[5] /. {SMP["d_g"] -> SMP["d_g^MS"]}
solMSbar4 = sol[5] /. {SMP["d_g"] -> SMP["d_g^MSbar"], 
       1/Epsilon -> SMP["Delta"]}
```

![1fvjxabo9cu7q](img/1fvjxabo9cu7q.svg)

![0cck89zt6dmsq](img/0cck89zt6dmsq.svg)

![0fgc2k4cbtb03](img/0fgc2k4cbtb03.svg)

## Check the final results

```mathematica
knownResult = Factor2[{SMP["d_psi^MS"] -> 
           (-SMP["alpha_s"]/(4*Pi))*(1/Epsilon)*CF*GaugeXi["G"], 
         SMP["d_m^MS"] -> (-SMP["alpha_s"]/(4*Pi))*(1/Epsilon)*3*
             CF, SMP["d_psi^MSbar"] -> (-SMP["alpha_s"]/(4*Pi))*
             SMP["Delta"]*CF*GaugeXi["G"], SMP["d_m^MSbar"] -> 
           (-SMP["alpha_s"]/(4*Pi))*SMP["Delta"]*3*CF, 
         SMP["d_A^MS"] -> (SMP["alpha_s"]/(4*Pi))*(1/Epsilon)*
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
     knownResult, Text -> {"\tCompare to Muta, Foundations of \
    QCD, Eqs 2.5.131-2.5.147:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![04djx03zyqbxp](img/04djx03zyqbxp.svg)

![1gvc45mqkq8hn](img/1gvc45mqkq8hn.svg)