---
title: Neutrino down-type quark annihilation into a lepton and an up-type quark
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Nle Qdt -> Le Qut, EW, total cross section, tree"; 
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
MakeBoxes[pu, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(u\)]\)"; 
MakeBoxes[pd, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(d\)]\)"; 
MakeBoxes[pn, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(n\)]\)"; 
MakeBoxes[pl, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(l\)]\)"; 
```

Enable CKM mixing

```mathematica
$CKM = True; 
```

To avoid dealing with Goldstone bosons we do  the computation in the unitary gauge

```mathematica
InitializeModel[{SM, UnitarySM}, GenericModel -> 
       {Lorentz, UnitaryLorentz}]; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {F[4, {1, 1}], F[1, {2}]} -> {F[3, {1, 1}], F[2, {2}]}, 
       InsertionLevel -> {Classes}, Model -> {SM, UnitarySM}, 
       GenericModel -> {Lorentz, UnitaryLorentz}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0l3n3g08hb0aw](img/0l3n3g08hb0aw.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, 
       GaugeRules -> {FAGaugeXi[W | Z] -> Infinity}], 
     IncomingMomenta -> {pd, pn}, OutgoingMomenta -> {pu, pl}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![076ju0oxic6us](img/076ju0oxic6us.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, pd, pn, -pu, -pl, SMP["m_d"], 0, 
     SMP["m_u"], SMP["m_l"]]; 
```

## Square the amplitude

There is no polarization averaging for neutrinos here, as right handed neutrinos do not interact

```mathematica
ampSquared[0] = Factor[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
         amp[0]*ComplexConjugate[amp[0]]]]]
```

![0cd5e83gght31](img/0cd5e83gght31.svg)

In the following we neglect the momentum in the W-propagator as compared to the W-mass. This is a very good approximation at low energies, as then (pl-pn)^2  <= m_mu^2 << m_W^2.

```mathematica
ampSquared[1] = Normal[
     (Series[#1, {SMP["m_W"], Infinity, 0}] & )[
       FeynAmpDenominatorExplicit[(#1 /. {pl - pn -> 0} & )[
           FCE[ampSquared[0]]]]]]
```

![183i59xppao0x](img/183i59xppao0x.svg)

## Total cross section

The total cross-section 

```mathematica
prefac = 4*(Pi/(64*Pi^2*s))*
     (Sqrt[(s - SMP["m_l"]^2 - SMP["m_u"]^2)^2 - 
            4*SMP["m_l"]^2*SMP["m_u"]^2]/Sqrt[(s - SMP["m_d"]^2)^2])
```

![1t0o2nmmmol46](img/1t0o2nmmmol46.svg)

```mathematica
crossSectionTotal = PowerExpand[prefac*ampSquared[1]]
```

![0gf0rbhogb9ix](img/0gf0rbhogb9ix.svg)

If the lepton is a muon, then the up quark mass can be neglected

```mathematica
crossSectionTotalMuon = PowerExpand[crossSectionTotal /. 
       {SMP["m_u"] -> 0, SMP["m_l"] -> SMP["m_mu"]}]
```

![01rbj2w91t9qz](img/01rbj2w91t9qz.svg)

## Check the final results

```mathematica
knownResults
```

![1w1ptnvz2hyxg](img/1w1ptnvz2hyxg.svg)

```mathematica
knownResults = {(SMP["G_F"]^2*(s - SMP["m_mu"]^2)^2*
            SMP["V_ud", -I]*SMP["V_ud", I])/(Pi*s)}; 
FCCompareResults[{crossSectionTotalMuon}, knownResults, 
     Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
    Physics, Chapter 5.3:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![0mmfttuonqlxc](img/0mmfttuonqlxc.svg)

![1sn90nppe3tk8](img/1sn90nppe3tk8.svg)