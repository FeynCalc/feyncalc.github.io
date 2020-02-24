---
title: W decaying into an electron and an electron-antineutrino
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "W -> El Anel, EW, total decay rate, tree"; 
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

Nicer typesetting

```mathematica
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[3]} -> {F[2, {1}], -F[1, {1}]}, InsertionLevel -> 
         {Particles}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0g3dgblqyweg7](img/0g3dgblqyweg7.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![05otfzq8z8yf2](img/05otfzq8z8yf2.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]
SP[p] = SMP["m_W"]^2; 
SP[k1] = SMP["m_e"]^2; 
SP[k2] = 0; 
SP[k1, k2] = (SMP["m_W"]^2 - SMP["m_e"]^2)/2; 
SP[p, k1] = (SMP["m_W"]^2 + SMP["m_e"]^2)/2; 
SP[p, k2] = (SMP["m_W"]^2 - SMP["m_e"]^2)/2; 
```

## Square the amplitude

We average over the polarizations of the W-boson, hence the additional factor 1/3

```mathematica
ampSquared[0] = Simplify[
     (DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       DiracSimplify[FermionSpinSum[SUNSimplify[
             amp[0]*ComplexConjugate[amp[0]]]]]]]
```

![01bxlau5ouk48](img/01bxlau5ouk48.svg)

## Total decay rate

```mathematica
phaseSpacePrefactor[m1_, m2_, M_] := (1/(16*Pi*M))*
       Sqrt[1 - (m1 + m2)^2/M^2]*Sqrt[1 - (m1 - m2)^2/M^2]; 
```

```mathematica
totalDecayRate = Simplify[phaseSpacePrefactor[SMP["m_e"], 0, 
         SMP["m_W"]]*ampSquared[0]]
```

![0djv5seudv6e8](img/0djv5seudv6e8.svg)

## Check the final results

```mathematica
knownResults = {(SMP["G_F"]*(SMP["m_e"] - SMP["m_W"])^2*
            (SMP["m_e"] + SMP["m_W"])^2*(SMP["m_e"]^2 + 
               2*SMP["m_W"]^2))/(12*Sqrt[2]*Pi*SMP["m_W"]^3)}; 
FCCompareResults[{totalDecayRate}, knownResults, 
     Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
    Physics, Chapter 5.2:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1aq8ptrqt9wgt](img/1aq8ptrqt9wgt.svg)

![0q5so78n5fptx](img/0q5so78n5fptx.svg)