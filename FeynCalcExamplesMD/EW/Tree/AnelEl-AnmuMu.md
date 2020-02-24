---
title: Electron electron-antineutrino annihilation into a muon-antineutrino and a muon
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Anel El -> Anmu Mu, EW, total cross section, tree"; 
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
MakeBoxes[pe, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(e\)]\)"; 
MakeBoxes[pm, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(m\)]\)"; 
```

To avoid dealing with Goldstone bosons we do  the computation in the unitary gauge

```mathematica
InitializeModel[{SM, UnitarySM}, GenericModel -> 
       {Lorentz, UnitaryLorentz}]; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {-F[1, {1}], F[2, {1}]} -> {-F[1, {2}], F[2, {2}]}, 
       InsertionLevel -> {Classes}, Model -> {SM, UnitarySM}, 
       GenericModel -> {Lorentz, UnitaryLorentz}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![17d6oqvsala3z](img/17d6oqvsala3z.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, 
       GaugeRules -> {FAGaugeXi[W | Z] -> Infinity}], 
     IncomingMomenta -> {q1, pe}, OutgoingMomenta -> {q2, pm}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![0u0xxu0c33iyt](img/0u0xxu0c33iyt.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, q1, pe, -q2, -pm, 0, SMP["m_e"], 0, 
     SMP["m_mu"]]; 
```

## Square the amplitude

There is no polarization averaging for neutrinos here, as right handed neutrinos do not interact

```mathematica
ampSquared[0] = Factor[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
         amp[0]*ComplexConjugate[amp[0]]]]]
```

![1dc483s0mmpja](img/1dc483s0mmpja.svg)

In the following we neglect the momentum in the W-propagator as compared to the W-mass. This is a very good approximation at low energies, as then (pm-q2)^2  <= m_mu^2 << m_W^2.

```mathematica
ampSquared[1] = Normal[
     (Series[#1, {SMP["m_W"], Infinity, 0}] & )[
       FeynAmpDenominatorExplicit[(#1 /. {pm + q2 -> 0} & )[
           FCE[ampSquared[0]]]]]]
```

![1ew8nezilsozj](img/1ew8nezilsozj.svg)

## Total cross section

We need to carry out the angular integration, so let us specify the values of the temporal and spatial components of the 4-vectors

```mathematica
TC[q2] = (s - SMP["m_mu"]^2)/(2*Sqrt[s]); 
TC[pe] = (s + SMP["m_e"]^2)/(2*Sqrt[s]); 
CSP[pe] = (s - SMP["m_e"]^2)^2/(4*s); 
CSP[pm] = (s - SMP["m_mu"]^2)^2/(4*s); 
CSP[q2] = (s - SMP["m_mu"]^2)^2/(4*s); 
```

```mathematica
prefac = 2*(Pi/(64*Pi^2*s))*(Sqrt[(s - SMP["m_mu"]^2)^2]/
        Sqrt[(s - SMP["m_e"]^2)^2])
```

![1bxyp3vw6sj7o](img/1bxyp3vw6sj7o.svg)

```mathematica
integral = Integrate[Simplify[ampSquared[1] /. 
         u -> SMP["m_e"]^2 - 2*(TC[q2]*TC[pe] - Sqrt[CSP[q2]]*
                    Sqrt[CSP[pe]]*x)], {x, -1, 1}]
```

![0c6jxr8g9dlut](img/0c6jxr8g9dlut.svg)

The total cross-section 

```mathematica
crossSectionTotal = Factor2[PowerExpand[integral*prefac]]
```

![151w72q75592q](img/151w72q75592q.svg)

## Check the final results

```mathematica
knownResults = {(SMP["G_F"]^2*(s - SMP["m_mu"]^2)^2)*
         ((s^2 + ((SMP["m_mu"]^2 + SMP["m_e"]^2)/2)*s + 
               SMP["m_mu"]^2*SMP["m_e"]^2)/(3*Pi*s^3))}; 
FCCompareResults[{crossSectionTotal}, knownResults, 
     Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
    Physics, Chapter 5.3:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![1j11dl44ry6v4](img/1j11dl44ry6v4.svg)

![0sqzwz0ujaboi](img/0sqzwz0ujaboi.svg)