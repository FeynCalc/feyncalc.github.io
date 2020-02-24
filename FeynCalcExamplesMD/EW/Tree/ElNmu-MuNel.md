---
title: Electron muon-neutrino annihilation into a muon and an electron-neutrino
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Nmu -> Mu Nuel, EW, total cross section, tree"; 
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
MakeBoxes[q1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(q\), \(1\)]\)"; 
MakeBoxes[q2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(q\), \(2\)]\)"; 
```

To avoid dealing with Goldstone bosons we do  the computation in the unitary gauge

```mathematica
InitializeModel[{SM, UnitarySM}, GenericModel -> 
       {Lorentz, UnitaryLorentz}]; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {F[2, {1}], F[1, {2}]} -> {F[1, {1}], F[2, {2}]}, 
       InsertionLevel -> {Classes}, Model -> {SM, UnitarySM}, 
       GenericModel -> {Lorentz, UnitaryLorentz}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1ln56t954awqf](img/1ln56t954awqf.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, 
       GaugeRules -> {FAGaugeXi[W | Z] -> Infinity}], 
     IncomingMomenta -> {p, q1}, OutgoingMomenta -> {q2, k}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![0q9hrkiuj6xce](img/0q9hrkiuj6xce.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p, q1, -q2, -k, SMP["m_e"], 0, 0, 
     SMP["m_mu"]]; 
```

## Square the amplitude

There is no polarization averaging for neutrinos here, as right handed neutrinos do not interact

```mathematica
ampSquared[0] = Factor[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
         amp[0]*ComplexConjugate[amp[0]]]]]
```

![0wc3du9yttix2](img/0wc3du9yttix2.svg)

In the following we neglect the momentum in the W-propagator as compared to the W-mass. This is a very good approximation at low energies, as then (k-q1)^2  <= m_mu^2 << m_W^2.

```mathematica
ampSquared[1] = Normal[
     (Series[#1, {SMP["m_W"], Infinity, 0}] & )[
       FeynAmpDenominatorExplicit[(#1 /. {k - q1 -> 0} & )[
           FCE[ampSquared[0]]]]]]
```

![1kcf6latg8he3](img/1kcf6latg8he3.svg)

## Total cross section

The total cross-section 

```mathematica
prefac = 4*(Pi/(64*Pi^2*s))*(Sqrt[(s - SMP["m_mu"]^2)^2]/
        Sqrt[(s - SMP["m_e"]^2)^2])
```

![0a2qe9fq3kjib](img/0a2qe9fq3kjib.svg)

```mathematica
crossSectionTotal = PowerExpand[prefac*ampSquared[1]]
```

![1bgoo5uu9e5kq](img/1bgoo5uu9e5kq.svg)

## Check the final results

```mathematica
knownResults = {(SMP["G_F"]^2*(s - SMP["m_mu"]^2)^2)/(Pi*s)}; 
FCCompareResults[{crossSectionTotal}, knownResults, 
     Text -> {"\tCompare to Greiner and Mueller, Gauge Theory of \
    Weak Interactions, Chapter 3:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1azloqyo538vk](img/1azloqyo538vk.svg)

![03trej0mp7a4e](img/03trej0mp7a4e.svg)