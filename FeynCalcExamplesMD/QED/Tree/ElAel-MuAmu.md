---
title: Muon production
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Ael -> Mu Amu, QED, total cross section, tree"; 
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
MakeBoxes[p1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(1\)]\)"; 
MakeBoxes[p2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(2\)]\)"; 
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {F[2, {1}], -F[2, {1}]} -> {F[2, {2}], -F[2, {2}]}, 
       InsertionLevel -> {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![12p4t3pmnxuzy](img/12p4t3pmnxuzy.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True]
```

![06vfbebevqwr5](img/06vfbebevqwr5.svg)

Polarized production: the particles are right-handed, while the antiparticles are left-handed

```mathematica
ampPolarized[0] = amp[0] /. {Spinor[-Momentum[k2], r__] :> 
         GA[6] . Spinor[-Momentum[k2], r], 
       Spinor[Momentum[k1], r__] :> Spinor[Momentum[k1], r] . 
           GA[7], Spinor[Momentum[p1], r__] :> 
         GA[6] . Spinor[Momentum[p1], r], 
       Spinor[-Momentum[p2], r__] :> Spinor[Momentum[p2], r] . 
           GA[7]}
```

![1lhyl87ljv529](img/1lhyl87ljv529.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"], 
     SMP["m_e"], SMP["m_mu"], SMP["m_mu"]]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[DiracSimplify[
       (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
         FeynAmpDenominatorExplicit[amp[0]*ComplexConjugate[
               amp[0]]]]]]
```

![0tvdxr7nkubvx](img/0tvdxr7nkubvx.svg)

```mathematica
ampSquaredPolarized[0] = Simplify[DiracSimplify[
       FermionSpinSum[FeynAmpDenominatorExplicit[
           ampPolarized[0]*ComplexConjugate[ampPolarized[0]]]]]]
```

![06apy7izzj4jl](img/06apy7izzj4jl.svg)

```mathematica
ampSquaredMassless1[0] = (#1 /. {SMP["m_e"] -> 0} & )[
     ampSquared[0]]
```

![0u7nu5yf296eu](img/0u7nu5yf296eu.svg)

```mathematica
ampSquaredMassless2[0] = Simplify[
     (#1 /. {SMP["m_e"] -> 0, SMP["m_mu"] -> 0} & )[
       ampSquared[0]]]
```

![1h25tuebaxog2](img/1h25tuebaxog2.svg)

```mathematica
ampSquaredPolarizedMassless[0] = 
   Simplify[(#1 /. {SMP["m_e"] -> 0, SMP["m_mu"] -> 0} & )[
       ampSquaredPolarized[0]]]
```

![0lzbnjcaat2u1](img/0lzbnjcaat2u1.svg)

## Total cross-section

The differential cross-section d sigma/ d Omega is given by

```mathematica
prefac1 = 1/(64*Pi^2*s); 
```

```mathematica
integral1 = Factor[ampSquaredMassless2[0] /. 
       {t -> (-s/2)*(1 - Cos[Th]), u -> (-s/2)*(1 + Cos[Th]), 
         SMP["e"]^4 -> (4*Pi*SMP["alpha_fs"])^2}]
```

![1wsvq0tmw9ayd](img/1wsvq0tmw9ayd.svg)

```mathematica
diffXSection1 = prefac1*integral1
```

![1bs8e59snlon6](img/1bs8e59snlon6.svg)

The differential cross-section d sigma/ d t d phi is given by

```mathematica
prefac2 = 1/(128*Pi^2*s)
```

![1iod6wvzt9h6n](img/1iod6wvzt9h6n.svg)

```mathematica
integral2 = Simplify[ampSquaredMassless2[0]/(s/4) /. 
       {u -> -s - t, SMP["e"]^4 -> (4*Pi*SMP["alpha_fs"])^2}]
```

![1thp1llce562f](img/1thp1llce562f.svg)

```mathematica
diffXSection2 = prefac2*integral2
```

![0bcsr5mrz979o](img/0bcsr5mrz979o.svg)

The total cross-section. We see that integrating both expressions gives the same result

```mathematica
2*Pi*Integrate[diffXSection1*Sin[Th], {Th, 0, Pi}]
```

![0i7tltqrnen5b](img/0i7tltqrnen5b.svg)

```mathematica
crossSectionTotal = 2*Pi*Integrate[diffXSection2, {t, -s, 0}]
```

![1d0c90ibcv1lv](img/1d0c90ibcv1lv.svg)

## Check the final results

```mathematica
knownResults = {(#1 /. SMP["m_e"] -> 0 & )[
         ExpandScalarProduct[(8*SMP["e"]^4*(SP[p1, k1]*SP[p2, k2] + 
                   SP[p1, k2]*SP[p2, k1] + 
          SMP["m_mu"]^2*SP[p1, p2]))/
             SP[p1 + p2]^2]], (16*
      SMP["e"]^4*(SP[p1, k2]*SP[p2, k1]))/
         SP[p1 + p2]^2, (8*(SMP["e"]^4/s^2))*((t/2)^2 + (u/2)^2), 
       (4*Pi*SMP["alpha_fs"]^2)/(3*s)}; 
FCCompareResults[{ampSquaredMassless1[0], 
       ampSquaredPolarized[0], ampSquaredMassless2[0], 
       crossSectionTotal}, knownResults, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eqs 5.10, 5.21, 5.70 and to Field, Applications of \
    Perturbative QCD, Eq. 2.1.14", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0w6rxpg8imyu8](img/0w6rxpg8imyu8.svg)

![05m28k1f06g9b](img/05m28k1f06g9b.svg)