---
title: Electron electron-antineutrino annihilation into an antiup and a down quark
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Anel El -> Qubar Qd, EW, total cross section, tree"; 
If[$FrontEnd === Null, $FeynCalcStartupMessages = False; 
      Print[description]; ]; 
If[$Notebooks === False, $FeynCalcStartupMessages = False]; 
$LoadAddOns = {"FeynArts"}; 
Get["FeynCalc`"]
$FAVerbose = 0; 
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
       {-F[1, {1}], F[2, {1}]} -> {-F[3, {1}], F[4, {1}]}, 
       InsertionLevel -> {Classes}, Model -> {SM, UnitarySM}, 
       GenericModel -> {Lorentz, UnitaryLorentz}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1ssn0hhn3fi90](img/1ssn0hhn3fi90.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, 
       GaugeRules -> {FAGaugeXi[W | Z] -> Infinity}], 
     DropSumOver -> True, IncomingMomenta -> {q, pe}, 
     OutgoingMomenta -> {l1, l2}, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![04xxzbwqejl2e](img/04xxzbwqejl2e.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, q, pe, -l1, -l2, 0, SMP["m_e"], 
     SMP["m_u"], SMP["m_d"]]; 
```

## Square the amplitude

There is no polarization averaging for neutrinos here, as right handed neutrinos do not interact

```mathematica
ampSquared[0] = Factor[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
         amp[0]*ComplexConjugate[amp[0]]]]]
```

![1h7ohpplu7t43](img/1h7ohpplu7t43.svg)

In the following we neglect the momentum in the W-propagator as compared to the W-mass. This is a very good approximation at low energies, as then (pm-q2)^2  <= m_mu^2 << m_W^2.

```mathematica
ampSquared[1] = Normal[
     (Series[#1, {SMP["m_W"], Infinity, 0}] & )[
       FeynAmpDenominatorExplicit[(#1 /. {l1 + l2 -> 0} & )[
           FCE[ampSquared[0]]]]]]
```

![0quvrh08b3wvt](img/0quvrh08b3wvt.svg)

## Total cross-section

We need to carry out the angular integration, so let us specify the values of the temporal and spatial components of the 4-vectors

```mathematica
TC[pe] = (s + SMP["m_e"]^2)/(2*Sqrt[s]); 
TC[l1] = (s + (SMP["m_u"]^2 - SMP["m_d"]^2))/(2*Sqrt[s]); 
CSP[pe] = (s - SMP["m_e"]^2)^2/(4*s); 
CSP[l1] = ((s - SMP["m_u"]^2 - SMP["m_d"]^2)^2 - 
          4*SMP["m_u"]^2*SMP["m_d"]^2)/(4*s); 
```

```mathematica
prefac = 2*(Pi/(64*Pi^2*s))*(Sqrt[CSP[l1]]/Sqrt[CSP[pe]]); 
```

```mathematica
integral = Simplify[SUNSimplify[
       FRH[(Integrate[#1, {x, -1, 1}] & )[
           Collect2[ampSquared[1] /. u -> SMP["m_e"]^2 + 
                   SMP["m_u"]^2 - 2*(TC[l1]*TC[pe] - Sqrt[CSP[l1]]*
                          Sqrt[CSP[pe]]*x), x, IsolateNames -> KK]]]]]
```

![1c7wxk7742ait](img/1c7wxk7742ait.svg)

The total cross-section 

```mathematica
crossSectionTotal = Factor2[PowerExpand[integral*prefac]]
```

![067k67376k959](img/067k67376k959.svg)

## Check the final results

```mathematica
knownResults = 
     {(CA*SMP["G_F"]^2*Sqrt[(s - SMP["m_d"]^2 - 
                   2*SMP["m_d"]*SMP["m_u"] - SMP["m_u"]^2)*
                (s - SMP["m_d"]^2 + 2*SMP["m_d"]*SMP["m_u"] - 
                   SMP["m_u"]^2)]*(2*s^3 - s^2*SMP["m_d"]^2 - 
               s*SMP["m_d"]^4 + s^2*SMP["m_e"]^2 + s*SMP["m_d"]^2*
                 SMP["m_e"]^2 - 2*SMP["m_d"]^4*SMP["m_e"]^2 - 
               s^2*SMP["m_u"]^2 + 2*s*SMP["m_d"]^2*SMP["m_u"]^2 + 
               s*SMP["m_e"]^2*SMP["m_u"]^2 + 4*SMP["m_d"]^2*
                 SMP["m_e"]^2*SMP["m_u"]^2 - s*SMP["m_u"]^4 - 
               2*SMP["m_e"]^2*SMP["m_u"]^4)*SMP["V_ud", -I]*
            SMP["V_ud", I])/(6*Pi*s^3)}; 
FCCompareResults[{crossSectionTotal}, knownResults, 
     Text -> {"\tCompare to the known result:", "CORRECT.", 
         "WRONG!"}, Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![14u442h82bzym](img/14u442h82bzym.svg)

![1bc1gs268rgif](img/1bc1gs268rgif.svg)
