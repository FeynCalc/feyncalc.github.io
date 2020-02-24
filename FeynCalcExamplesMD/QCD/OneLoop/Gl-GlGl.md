---
title: QCD 3-gluon vertex at 1-loop
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom QCD model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/QCD/GenerateModelQCD.m before running it for the first time.

```mathematica
description = "Gl - Gl Gl, QCD, only UV divergences, 1-loop"; 
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
```

```mathematica
template = insertFields[createTopologies[1, 1 -> 2, 
         ExcludeTopologies -> {Tadpoles, WFCorrections, 
             WFCorrectionCTs, SelfEnergies}], {V[5]} -> {V[5], 
     V[5]}, 
       InsertionLevel -> {Particles}, 
       Model -> FileNameJoin[{"QCD", "QCD"}], 
       GenericModel -> FileNameJoin[{"QCD", "QCD"}], 
       ExcludeParticles -> {F[3 | 4, {2 | 3}], F[4, {1}]}]; 
```

```mathematica
diags = template /. createTopologies -> CreateTopologies /. 
       insertFields -> InsertFields; 
diagsCT = template /. createTopologies -> CreateCTTopologies /. 
       insertFields -> InsertFields; 
```

```mathematica
Paint[diags, ColumnsXRows -> {3, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![14lp9mg4prcy3](img/14lp9mg4prcy3.svg)

![0arc1gna7gjxf](img/0arc1gna7gjxf.svg)

![143ujnsb3fabs](img/143ujnsb3fabs.svg)

```mathematica
Paint[diagsCT, ColumnsXRows -> {3, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![197ji0fw1gm26](img/197ji0fw1gm26.svg)

## Obtain the amplitudes

The 1/(2Pi)^D prefactor is implicit. We keep the full gauge dependence. To simplify comparisons
to the literature, we make all momenta incoming.

Quark contribution. Notice that we multiply the amplitude by Nf to account for the number
of quark flavours in the loop.

```mathematica
amp1[0] = FCFAConvert[CreateFeynAmp[DiagramExtract[diags, 
         {1, 2}], Truncated -> True, GaugeRules -> {}, 
       PreFactor -> 1], IncomingMomenta -> {p1}, 
     OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
       {mu, nu, rho}, DropSumOver -> True, 
     SUNIndexNames -> {a, b, c}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, Prefactor -> Nf, 
     FinalSubstitutions -> {SMP["m_u"] -> SMP["m_q"], p2 -> -p2, 
         p3 -> -p3}]
```

![1kfbjbn3qyc8h](img/1kfbjbn3qyc8h.svg)

Ghost contribution

```mathematica
amp2[0] = FCFAConvert[CreateFeynAmp[DiagramExtract[diags, 
         {3, 4}], Truncated -> True, GaugeRules -> {}, 
       PreFactor -> 1], IncomingMomenta -> {p1}, 
     OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
       {mu, nu, rho}, SUNIndexNames -> {a, b, c}, 
     DropSumOver -> True, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {SMP["m_u"] -> SMP["m_q"], p2 -> -p2, p3 -> -p3}]
```

![1tylvyuzlv2bw](img/1tylvyuzlv2bw.svg)

Gluon contribution

```mathematica
amp3[0] = FCFAConvert[CreateFeynAmp[DiagramExtract[diags, 
           {5, 6, 7, 8}], Truncated -> True, GaugeRules -> {}, 
         PreFactor -> 1], IncomingMomenta -> {p1}, 
       OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
         {mu, nu, rho}, SUNIndexNames -> {a, b, c}, 
       DropSumOver -> True, LoopMomenta -> {l}, 
       UndoChiralSplittings -> True, ChangeDimension -> D, 
       List -> False, SMP -> True, FinalSubstitutions -> 
         {SMP["m_u"] -> SMP["m_q"], p2 -> -p2, p3 -> -p3}]; 
```

Counter-term

```mathematica
amp4[0] = FCFAConvert[CreateFeynAmp[diagsCT, Truncated -> True, 
       GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p1}, 
     OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
       {mu, nu, rho}, SUNIndexNames -> {a, b, c}, 
     DropSumOver -> True, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {SMP["m_u"] -> SMP["m_q"], p2 -> -p2, p3 -> -p3, 
         ZA -> SMP["Z_A"], Zg -> SMP["Z_g"]}]
```

![1f67kjqpygr67](img/1f67kjqpygr67.svg)

## Calculate the amplitudes

### Quark contribution

```mathematica
AbsoluteTiming[
   amp1[1] = TID[FCE[amp1[0]] /. {p2 + p3 -> -p1, 
             -p2 - p3 -> p1}, l, UsePaVeBasis -> True, 
         ToPaVe -> True]; ]
```

![06yuq7f9w2k6b](img/06yuq7f9w2k6b.svg)

```mathematica
amp1Div[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D, 
            FCLoopExtract -> False] & )[amp1[1]]; 
```

```mathematica
amp1Div[1] = (Collect2[#1, MTD, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[(SelectNotFree2[#1, SMP["Delta"]] & )[
         FCHideEpsilon[Normal[(Series[#1, {Epsilon, 0, 0}] & )[
               (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
                 (#1 /. SUNTrace[x__] :> SUNTrace[x, Explicit -> 
                              True] & )[(SUNSimplify[#1, 
             Explicit -> True] & )[
                     amp1Div[0]]]]]]]]]]
```

![1oo12t7403hw3](img/1oo12t7403hw3.svg)

In the calculation p3 was eliminated via the 4-momentum conservation p1+p2+p3=0. 
Now we need to reintroduce it

```mathematica
amp1Div[2] = amp1Div[1] /. {2*p1 + p2 -> p1 - p3, 
       p1 + 2*p2 -> p2 - p3}
```

![01s8xgyo78h3q](img/01s8xgyo78h3q.svg)

### Ghost contribution

```mathematica
AbsoluteTiming[
   amp2[1] = TID[FCE[amp2[0]] /. {p2 + p3 -> -p1, 
             -p2 - p3 -> p1}, l, UsePaVeBasis -> True, 
         ToPaVe -> True]; ]
```

![0utodty6w1dh5](img/0utodty6w1dh5.svg)

```mathematica
amp2Div[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D, 
            FCLoopExtract -> False] & )[amp2[1]]; 
```

```mathematica
amp2Div[1] = (Collect2[#1, MTD, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[(SelectNotFree2[#1, SMP["Delta"]] & )[
         FCHideEpsilon[Normal[(Series[#1, {Epsilon, 0, 0}] & )[
               (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
                 (#1 /. SUNTrace[x__] :> SUNTrace[x, Explicit -> 
                              True] & )[(SUNSimplify[#1, 
             Explicit -> True] & )[
                     amp2Div[0]]]]]]]]]]
```

![0x5lqlekfzapg](img/0x5lqlekfzapg.svg)

In the calculation p3 was eliminated via the 4-momentum conservation p1+p2+p3=0. 
Now we need to reintroduce it

```mathematica
amp2Div[2] = amp2Div[1] /. {2*p1 + p2 -> p1 - p3, 
       p1 + 2*p2 -> p2 - p3}
```

![0vse4837ikswl](img/0vse4837ikswl.svg)

### Gluon contribution

This calculation requires about 70 seconds on a modern laptop

```mathematica
AbsoluteTiming[
   amp3[1] = TID[FCE[amp3[0]] /. {p2 + p3 -> -p1, 
             -p2 - p3 -> p1}, l, UsePaVeBasis -> True, 
         ExpandScalarProduct -> False, ToPaVe -> True]; ]
```

![1eh2wnrw0rj0g](img/1eh2wnrw0rj0g.svg)

```mathematica
amp3Div[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D, 
            FCLoopExtract -> False] & )[amp3[1]]; 
```

```mathematica
amp3Div[1] = (Collect2[#1, MTD, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[(SelectNotFree2[#1, SMP["Delta"]] & )[
         FCHideEpsilon[Normal[(Series[#1, {Epsilon, 0, 0}] & )[
               (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[
                 SUNSimplify[amp3Div[0]]]]]]]]]
```

![1a3w6rgl45nv2](img/1a3w6rgl45nv2.svg)

```mathematica
amp3Div[2] = (Collect2[#1, MTD, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[ExpandScalarProduct[
         (Collect2[#1, MTD, GaugeXi, Factoring -> Function[x, 
                            MomentumCombine[Factor[x]]]] & )[
                   FCE[ExpandScalarProduct[amp3Div[1] /. 
                         {p3 -> -p1 - p2}]]] /. 
        2*p1 + p2 -> p1 - p3 /. 
               p1 + 2*p2 -> p2 - p3 /. 22*p1 + 11*p2 -> 
               11*p1 - 11*p3 /. -22*p2 - 11*p1 -> -11*p2 + 11*p3]]]
```

![12pwnfpyngpgi](img/12pwnfpyngpgi.svg)

### Counter-term

```mathematica
amp4[1] = (Collect2[#1, MTD, GaugeXi, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[ExpandScalarProduct[(#1 /. alpha -> 1 & )[
           Normal[(Series[#1, {alpha, 0, 1}] & )[
               (#1 /. {SMP["Z_A"] -> 1 + alpha*SMP["d_A"], 
                        SMP["Z_g"] -> 1 + alpha*SMP["d_g"]} & )[
                 amp4[0]]]]]]]]
```

![0lrhkqv0de7dd](img/0lrhkqv0de7dd.svg)

Check the cancellation of the UV divergences in the MSbar scheme. The renormalization constants
are obtained from another example calculation, "Renormalization.m"

```mathematica
renormalizationConstants = 
     {SMP["d_A"] -> (SMP["alpha_s"]/(4*Pi))*SMP["Delta"]*
             ((1/2)*CA*(13/3 - GaugeXi["G"]) - (2/3)*Nf), 
         SMP["d_g"] -> ((-11*CA*SMP["alpha_s"])/(24*Pi))*
               SMP["Delta"] + ((Nf*SMP["alpha_s"])/(12*Pi))*
               SMP["Delta"]} /. 
   SMP["alpha_s"] -> SMP["g_s"]^2/(4*Pi); 
```

```mathematica
uvDiv[0] = (Collect2[#1, MTD, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[ExpandScalarProduct[amp1Div[2] + amp2Div[2] + 
           amp3Div[2] + amp4[1]]]]
```

![0dgrr1lhv2mez](img/0dgrr1lhv2mez.svg)

```mathematica
uvDiv[1] = Simplify[uvDiv[0] /. renormalizationConstants]
```

![1qih75y48ehqj](img/1qih75y48ehqj.svg)

```mathematica
FCCompareResults[uvDiv[1], 0, Text -> {"\tThe UV divergence of \
    the 3-gluon vertex at 1-loop is cancelled by the counter-term \
    :", "CORRECT.", "WRONG!"}, Interrupt -> {Hold[Quit[1]], 
         Automatic}]; 
```

![0h6qxqe1xggwd](img/0h6qxqe1xggwd.svg)

## Check the final results

```mathematica
VertexLorentzStruct[{p_, q_, k_}, {mu_, nu_, si_}, 
       {a_, b_, c_}] := (-I)*SUNF[a, b, c]*
       (MTD[mu, nu]*FVD[p - q, si] + MTD[nu, si]*FVD[q - k, mu] + 
          MTD[si, mu]*FVD[k - p, nu]); 
```

```mathematica
knownResult = FCI[{(I*SMP["g_s"])*Nf*(SMP["g_s"]^2/(4*Pi)^2)*
           (-2/3)*SMP["Delta"]*VertexLorentzStruct[{p1, p2, p3}, 
             {mu, nu, rho}, {a, b, c}], (I*SMP["g_s"])*
           (SMP["g_s"]^2/(4*Pi)^2)*(CA/8)*(1/3)*SMP["Delta"]*
           VertexLorentzStruct[{p1, p2, p3}, {mu, nu, rho}, 
             {a, b, c}], (I*SMP["g_s"])*(SMP["g_s"]^2/(4*Pi)^2)*
           (CA/8)*SMP["Delta"]*(-4 - 9*GaugeXi["G"] + 15 + 
              3*GaugeXi["G"])*VertexLorentzStruct[{p1, p2, p3}, 
             {mu, nu, rho}, {a, b, c}]}]; 
FCCompareResults[{amp1Div[2], amp2Div[2], amp3Div[2]}, 
     knownResult, Text -> {"\tCompare to Pascual and Tarrach, \
    QCD: Renormalization for the Practitioner, Eq III.46:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![00bclrcuvh0hi](img/00bclrcuvh0hi.svg)

![03i9ef9w01j6j](img/03i9ef9w01j6j.svg)