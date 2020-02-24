---
title: Higgs decaying into a W boson pair
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "H -> W W, EW, total decay rate, tree"; 
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
       {S[1]} -> {-V[3], V[3]}, InsertionLevel -> {Classes}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0y6mt6ddwxfhx](img/0y6mt6ddwxfhx.svg)

## Obtain the amplitudes

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {pH}, OutgoingMomenta -> {k1, k2}, 
     List -> False, ChangeDimension -> 4, DropSumOver -> True, 
     SMP -> True, Contract -> True, UndoChiralSplittings -> True]
```

![0wicd200nputn](img/0wicd200nputn.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[k1, k1] = SMP["m_W"]^2; 
SP[k2, k2] = SMP["m_W"]^2; 
SP[pH, pH] = SMP["m_H"]^2; 
SP[k1, k2] = (SMP["m_H"]^2 - 2*SMP["m_W"]^2)/2; 
```

## Square the amplitudes

```mathematica
ampSquared[0] = Simplify[(DoPolarizationSums[#1, k2] & )[
       (DoPolarizationSums[#1, k1] & )[FeynAmpDenominatorExplicit[
           amp[0]*ComplexConjugate[amp[0]]]]]]
```

![19f9t45fww1uz](img/19f9t45fww1uz.svg)

## Total decay rate

```mathematica
$Assumptions = {SMP["m_H"] > 0, SMP["m_W"] > 0}; 
phaseSpacePrefactor[m_] := (1/(16*Pi*SMP["m_H"]))*
       Sqrt[1 - 4*(m^2/SMP["m_H"]^2)]; 
```

```mathematica
totalDecayRate = Simplify[
     (#1 /. SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"] & )[
       phaseSpacePrefactor[SMP["m_W"]]*ampSquared[0]]]
```

![08n9g6hkuafs9](img/08n9g6hkuafs9.svg)

Rewrite the result in a nicer way

```mathematica
(#1 /. h -> Identity & )[FullSimplify[
     totalDecayRate /. SMP["m_W"]^2 -> 
           h[SMP["m_W"]^2/SMP["m_H"]^2]*SMP["m_H"]^2 /. 
       SMP["m_W"]^4 -> h[SMP["m_W"]^4/SMP["m_H"]^4]*SMP["m_H"]^4]]
```

![0mcyvyxqcvw3c](img/0mcyvyxqcvw3c.svg)

## Check the final results

```mathematica
knownResults = (#1 /. SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"] & )[
       (#1 /. SMP["g_W"] -> SMP["e"]/SMP["sin_W"] & )[
         {(SMP["g_W"]^2/(64*Pi))*(SMP["m_H"]^3/SMP["m_W"]^2)*
             Sqrt[1 - 4*(SMP["m_W"]^2/SMP["m_H"]^2)]*
             (1 - 4*(SMP["m_W"]^2/SMP["m_H"]^2) + 
                12*(SMP["m_W"]^4/SMP["m_H"]^4))}]]; 
FCCompareResults[{totalDecayRate}, knownResults, 
     Factoring -> Simplify, Text -> {"\tCompare to Gunion, Haber, \
    Kane and Dawson, Higgs Hunter Guide, Eq 2.11:", "CORRECT.", 
         "WRONG!"}, Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![0b0ph2ooqhrel](img/0b0ph2ooqhrel.svg)

![0qef73c152x8f](img/0qef73c152x8f.svg)