---
title: Higgs decaying into a fermion-antifermion pair
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "H -> F Fbar, EW, total decay rate, tree"; 
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

```mathematica
diagsLeptons = InsertFields[CreateTopologies[0, 1 -> 2], 
       {S[1]} -> {F[2, {l}], -F[2, {l}]}, InsertionLevel -> 
         {Classes}]; 
Paint[diagsLeptons, ColumnsXRows -> {2, 1}, 
     Numbering -> Simple, SheetHeader -> None, 
     ImageSize -> {512, 256}]; 
```

![0uogn2z1390pm](img/0uogn2z1390pm.svg)

```mathematica
diagsQuarks = InsertFields[CreateTopologies[0, 1 -> 2], 
       {S[1]} -> {F[3, {l}], -F[3, {l}]}, InsertionLevel -> 
         {Classes}]; 
Paint[diagsQuarks, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![17jaurmpdfxpr](img/17jaurmpdfxpr.svg)

## Obtain the amplitudes

```mathematica
ampLeptons[0] = FCFAConvert[CreateFeynAmp[diagsLeptons], 
     IncomingMomenta -> {pH}, OutgoingMomenta -> {p1, p2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, UndoChiralSplittings -> True, 
     FinalSubstitutions -> {MLE[l] -> SMP["m_l"]}]
```

![0wq1z6xfuy9fm](img/0wq1z6xfuy9fm.svg)

```mathematica
ampQuarks[0] = FCFAConvert[CreateFeynAmp[diagsQuarks], 
     IncomingMomenta -> {pH}, OutgoingMomenta -> {k1, k2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, UndoChiralSplittings -> True, 
     FinalSubstitutions -> {MQU[l] -> SMP["m_q"]}]
```

![15mvdqyw78epg](img/15mvdqyw78epg.svg)

```mathematica
ampLeptons[1] = DiracSimplify[ampLeptons[0]]
```

![1d2nfu7gqrb1l](img/1d2nfu7gqrb1l.svg)

```mathematica
ampQuarks[1] = DiracSimplify[ampQuarks[0]]
```

![1hw677xv7g8nb](img/1hw677xv7g8nb.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[p1, p1] = SMP["m_l"]^2; 
SP[k1, k1] = SMP["m_q"]^2; 
SP[p2, p2] = SMP["m_l"]^2; 
SP[k2, k2] = SMP["m_q"]^2; 
SP[pH, pH] = SMP["m_H"]^2; 
SP[p1, p2] = (SMP["m_H"]^2 - 2*SMP["m_l"]^2)/2; 
SP[k1, k2] = (SMP["m_H"]^2 - 2*SMP["m_q"]^2)/2; 
```

## Square the amplitudes

```mathematica
ampLeptonsSquared[0] = Simplify[DiracSimplify[
       FermionSpinSum[ampLeptons[1]*ComplexConjugate[
             ampLeptons[1]]]]]
```

![0rrrq632ftoxe](img/0rrrq632ftoxe.svg)

```mathematica
ampQuarksSquared[0] = SUNSimplify[DiracSimplify[
       FermionSpinSum[ampQuarks[1]*ComplexConjugate[
             ampQuarks[1]]]]]
```

![0sxqca5zbtjnu](img/0sxqca5zbtjnu.svg)

## Total decay rates

```mathematica
$Assumptions = {SMP["m_H"] > 0, SMP["m_l"] > 0}; 
phaseSpacePrefactor[m_] := (1/(16*Pi*SMP["m_H"]))*
       Sqrt[1 - 4*(m^2/SMP["m_H"]^2)]; 
```

```mathematica
totalDecayRateLeptons = Simplify[
     (#1 /. SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"] & )[
       phaseSpacePrefactor[SMP["m_l"]]*ampLeptonsSquared[0]]]
```

![0dvg6rxd7lyr1](img/0dvg6rxd7lyr1.svg)

```mathematica
totalDecayRateQuarks = Simplify[
     (#1 /. SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"] & )[
       phaseSpacePrefactor[SMP["m_q"]]*ampQuarksSquared[0]]]
```

![0qs5gqtruzirj](img/0qs5gqtruzirj.svg)

```mathematica
Factor[((AlphaFS*SMP["m_H"])/(8*SMP["sin_W"]^2))*
       ((SMP["m_l"]^2/SMP["m_W"]^2)*
          (1 - 4*(SMP["m_l"]^2/SMP["m_H"]^2))^(3/2)) - 
     totalDecayRateLeptons]
```

![0fp8cboltbq1k](img/0fp8cboltbq1k.svg)

## Check the final results

```mathematica
knownResults = {((AlphaFS*SMP["m_H"])/(8*SMP["sin_W"]^2))*
         ((SMP["m_l"]^2/SMP["m_W"]^2)*
            (1 - 4*(SMP["m_l"]^2/SMP["m_H"]^2))^(3/2)), 
       (CA*SMP["alpha_fs"]*SMP["m_H"]*SMP["m_q"]^2*
            (1 - (4*SMP["m_q"]^2)/SMP["m_H"]^2)^(3/2))/
         (8*SMP["m_W"]^2*SMP["sin_W"]^2)}; 
FCCompareResults[{totalDecayRateLeptons, totalDecayRateQuarks}, 
     knownResults, Factoring -> Simplify, 
     Text -> {"\tCompare to Peskin and Schroeder,An Introduction \
    to QFT, Final Project III, part (a):", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1jz7ndv9iqvt8](img/1jz7ndv9iqvt8.svg)

![17sz8jzzhaz1j](img/17sz8jzzhaz1j.svg)