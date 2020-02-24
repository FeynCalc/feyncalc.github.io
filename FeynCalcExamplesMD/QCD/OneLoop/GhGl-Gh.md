---
title: QCD ghost-gluon vertex at 1-loop
---


## Load FeynCalc and the necessary add-ons or other packages

This example uses a custom QCD model created with FeynRules. Please evaluate the file
FeynCalc/Examples/FeynRules/QCD/GenerateModelQCD.m before running it for the first time.

```mathematica
description = "GhGl - Gh, QCD, only UV divergences, 1-loop"; 
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
             WFCorrectionCTs, SelfEnergies}], {U[5]} -> {V[5], 
     U[5]}, 
       InsertionLevel -> {Particles}, 
       Model -> FileNameJoin[{"QCD", "QCD"}], 
       GenericModel -> FileNameJoin[{"QCD", "QCD"}]]; 
```

```mathematica
diags = template /. createTopologies -> CreateTopologies /. 
       insertFields -> InsertFields; 
diagsCT = template /. createTopologies -> CreateCTTopologies /. 
       insertFields -> InsertFields; 
```

```mathematica
Paint[diags, ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![1va61gy7705zb](img/1va61gy7705zb.svg)

```mathematica
Paint[diagsCT, ColumnsXRows -> {2, 1}, SheetHeader -> None, 
     Numbering -> Simple, ImageSize -> {512, 256}]; 
```

![0l2ainpke3ts3](img/0l2ainpke3ts3.svg)

## Obtain the amplitudes

The 1/(2Pi)^D prefactor is implicit. We keep the full gauge dependence.

```mathematica
amp1[0] = FCFAConvert[CreateFeynAmp[DiagramExtract[diags, 
         {1, 2}], Truncated -> True, GaugeRules -> {}, 
       PreFactor -> 1], IncomingMomenta -> {p1}, 
     OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
       {mu, nu, rho}, DropSumOver -> True, 
     SUNIndexNames -> {a, b, c}, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {SMP["m_u"] -> SMP["m_q"]}]
```

![0yijw7gghatc1](img/0yijw7gghatc1.svg)

Counter-term

```mathematica
amp2[0] = FCFAConvert[CreateFeynAmp[diagsCT, Truncated -> True, 
       GaugeRules -> {}, PreFactor -> 1], IncomingMomenta -> {p1}, 
     OutgoingMomenta -> {p2, p3}, LorentzIndexNames -> 
       {mu, nu, rho}, SUNIndexNames -> {a, b, c}, 
     DropSumOver -> True, LoopMomenta -> {l}, 
     UndoChiralSplittings -> True, ChangeDimension -> D, 
     List -> False, SMP -> True, FinalSubstitutions -> 
       {SMP["m_u"] -> SMP["m_q"], ZA -> SMP["Z_A"], 
         Zg -> SMP["Z_g"], Zu -> SMP["Z_u"]}]
```

![0lh0iayycfr58](img/0lh0iayycfr58.svg)

## Calculate the amplitudes

### Ghost-gluon vertex

```mathematica
AbsoluteTiming[amp1[1] = TID[FCE[amp1[0]] /. {-p2 - p3 -> -p1}, 
         l, UsePaVeBasis -> True, ToPaVe -> True]; ]
```

![1apbqa6mxmbkv](img/1apbqa6mxmbkv.svg)

```mathematica
amp1Div[0] = (PaVeUVPart[#1, Prefactor -> 1/(2*Pi)^D] & )[
       amp1[1]]; 
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

![0cc2lh0wbzz9e](img/0cc2lh0wbzz9e.svg)

### Counter-term

```mathematica
amp2[1] = (Collect2[#1, MTD, GaugeXi, Factoring -> 
            Function[x, MomentumCombine[Factor[x]]]] & )[
     FCE[ExpandScalarProduct[(#1 /. alpha -> 1 & )[
           Normal[(Series[#1, {alpha, 0, 1}] & )[
               (#1 /. {SMP["Z_A"] -> 1 + alpha*SMP["d_A"], 
                        SMP["Z_u"] -> 1 + alpha*SMP["d_u"], 
            SMP["Z_g"] -> 
                          1 + alpha*SMP["d_g"]} & )[amp2[0]]]]]]]]
```

![1cgtp2eddle3o](img/1cgtp2eddle3o.svg)

Check the cancellation of the UV divergences in the MSbar scheme. The renormalization constants
are obtained from another example calculation, "Renormalization.m"

```mathematica
renormalizationConstants = 
     {SMP["d_A"] -> (SMP["alpha_s"]/(4*Pi))*SMP["Delta"]*
             ((1/2)*CA*(13/3 - GaugeXi["G"]) - (2/3)*Nf), 
         SMP["d_g"] -> ((-11*CA*SMP["alpha_s"])/(24*Pi))*
               SMP["Delta"] + ((Nf*SMP["alpha_s"])/(12*Pi))*
               SMP["Delta"], 
    SMP["d_u"] -> (SMP["alpha_s"]/(4*Pi))*CA*
             SMP["Delta"]*((3 - GaugeXi["G"])/4)} /. 
       SMP["alpha_s"] -> SMP["g_s"]^2/(4*Pi); 
```

```mathematica
uvDiv[0] = Simplify[ExpandScalarProduct[amp1Div[1] + amp2[1]]]
```

![0izgzts61gd5c](img/0izgzts61gd5c.svg)

```mathematica
uvDiv[1] = Simplify[uvDiv[0] /. renormalizationConstants]
```

![1b3x5u5rst8v1](img/1b3x5u5rst8v1.svg)

```mathematica
FCCompareResults[uvDiv[1], 0, Text -> {"\tThe UV divergence of \
    the ghost-gluon vertex at 1-loop is cancelled by the \
    counter-term :", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
```

![0uz4mx0gd6cjp](img/0uz4mx0gd6cjp.svg)

## Check the final results

```mathematica
knownResult = (-I)*((-I)*SMP["g_s"]*FVD[p3, mu]*SUNF[a, b, c]*
          ((SMP["g_s"]^2/(4*Pi)^2)*CA*(GaugeXi["G"]/2)*
             SMP["Delta"])); 
FCCompareResults[amp1Div[1], knownResult, 
     Text -> 
       {"\tCompare to Muta, Foundations of QCD, Eq. 2.5.142:", 
         "CORRECT.", "WRONG!"}, Interrupt -> 
       {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0zfzzhhj0lv9s](img/0zfzzhhj0lv9s.svg)

![0wv01mdr84rfe](img/0wv01mdr84rfe.svg)