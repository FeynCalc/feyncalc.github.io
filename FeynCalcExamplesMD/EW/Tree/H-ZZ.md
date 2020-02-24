---
title: Higgs decaying into a Z boson pair
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "H -> Z Z, EW, total decay rate, tree"; 
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
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {S[1]} -> {V[2], V[2]}, InsertionLevel -> {Classes}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0bei6n4xd4tdf](img/0bei6n4xd4tdf.svg)

## Obtain the amplitudes

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {pH}, OutgoingMomenta -> {k1, k2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, UndoChiralSplittings -> True]
```

![1cco6ob4fqrz3](img/1cco6ob4fqrz3.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[k1, k1] = SMP["m_Z"]^2; 
SP[k2, k2] = SMP["m_Z"]^2; 
SP[pH, pH] = SMP["m_H"]^2; 
SP[k1, k2] = (SMP["m_H"]^2 - 2*SMP["m_Z"]^2)/2; 
```

## Square the amplitudes

```mathematica
ampSquared[0] = Simplify[(DoPolarizationSums[#1, k2] & )[
       (DoPolarizationSums[#1, k1] & )[FeynAmpDenominatorExplicit[
           (1/2)*(amp[0]*ComplexConjugate[amp[0]])]]]]
```

![0njv8msto60s6](img/0njv8msto60s6.svg)

## Total decay rates

```mathematica
$Assumptions = {SMP["m_H"] > 0, SMP["m_Z"] > 0}; 
phaseSpacePrefactor[m_] := (1/(16*Pi*SMP["m_H"]))*
       Sqrt[1 - 4*(m^2/SMP["m_H"]^2)]; 
```

```mathematica
totalDecayRate = Simplify[
     (#1 //. {SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"], 
              1/SMP["m_Z"]^4 -> SMP["cos_W"]^4/SMP["m_W"]^4} & )[
       phaseSpacePrefactor[SMP["m_Z"]]*ampSquared[0]]]
```

![0izpapg7ks8nd](img/0izpapg7ks8nd.svg)

Rewrite the result in a nicer way

```mathematica
(#1 /. h -> Identity & )[FullSimplify[
     totalDecayRate /. SMP["m_Z"]^2 -> 
           h[SMP["m_Z"]^2/SMP["m_H"]^2]*SMP["m_H"]^2 /. 
       SMP["m_Z"]^4 -> h[SMP["m_Z"]^4/SMP["m_H"]^4]*SMP["m_H"]^4]]
```

![1ch4iyq3apx1m](img/1ch4iyq3apx1m.svg)

## Check the final results

```mathematica
knownResults = {(SMP["alpha_fs"]*SMP["m_H"]^3*
            Sqrt[1 - (4*SMP["m_Z"]^2)/SMP["m_H"]^2]*
            (1 - (4*SMP["m_Z"]^2)/SMP["m_H"]^2 + (12*SMP["m_Z"]^4)/
                 SMP["m_H"]^4))/(32*SMP["m_W"]^2*SMP["sin_W"]^2)}; 
FCCompareResults[{totalDecayRate}, knownResults, 
     Factoring -> Simplify, Text -> {"\tCompare to Gunion, Haber, \
    Kane and Dawson, Higgs Hunter Guide, Eq 2.10:", "CORRECT.", 
         "WRONG!"}, Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1n6dio04v9lqj](img/1n6dio04v9lqj.svg)

![1ump1w36iwq2x](img/1ump1w36iwq2x.svg)