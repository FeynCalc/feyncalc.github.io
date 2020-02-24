---
title: Quark-antiquark production in electron-positron annihilation
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Ael -> Q Qbar, QCD, total cross section, tree"; 
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
       {F[2, {1}], -F[2, {1}]} -> {F[3, {1}], -F[3, {1}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD", 
       ExcludeParticles -> {S[_], V[2]}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0t5kwehtwvu7k](img/0t5kwehtwvu7k.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     DropSumOver -> True, Prefactor -> (3/2)*SMP["e_Q"], 
     FinalSubstitutions -> {SMP["m_u"] -> SMP["m_q"]}]
```

![1xcqqtv2g3g3a](img/1xcqqtv2g3g3a.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"]^2, 
     SMP["m_e"]^2, SMP["m_q"]^2, SMP["m_q"]^2]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_q"]^2 + 
                2*SMP["m_e"]^2}] & )[DiracSimplify[
         (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                    False] & )[FeynAmpDenominatorExplicit[
               amp[0]*ComplexConjugate[amp[0]]]]]]]]
```

![09418tdvjkujn](img/09418tdvjkujn.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_q" | "m_e"] -> 0} & )[ampSquared[0]]]
```

![0mlek221ua9tu](img/0mlek221ua9tu.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![11lsfmp61x82y](img/11lsfmp61x82y.svg)

## Total cross-section

The differential cross-section d sigma/ d Omega is given by

```mathematica
prefac1 = 1/(64*Pi^2*s); 
```

```mathematica
integral1 = Factor[ampSquaredMasslessSUNN3[0] /. 
       {t -> (-s/2)*(1 - Cos[Th]), u -> (-s/2)*(1 + Cos[Th]), 
         SMP["e"]^4 -> (4*Pi*SMP["alpha_fs"])^2}]
```

![1i59hbkqbs9j0](img/1i59hbkqbs9j0.svg)

```mathematica
diffXSection1 = prefac1*integral1
```

![0brma2hlz8iey](img/0brma2hlz8iey.svg)

The differential cross-section d sigma/ d t d phi is given by

```mathematica
prefac2 = 1/(128*Pi^2*s)
```

![1rl95vxjvkz62](img/1rl95vxjvkz62.svg)

```mathematica
integral2 = Simplify[ampSquaredMasslessSUNN3[0]/(s/4) /. 
       {u -> -s - t, SMP["e"]^4 -> (4*Pi*SMP["alpha_fs"])^2}]
```

![08ds4gk2ddtlw](img/08ds4gk2ddtlw.svg)

```mathematica
diffXSection2 = prefac2*integral2
```

![03i8cmjrjypu2](img/03i8cmjrjypu2.svg)

The total cross-section. We see that integrating both expressions gives the same result

```mathematica
2*Pi*Integrate[diffXSection1*Sin[Th], {Th, 0, Pi}]
```

![0j8yme3vszlrg](img/0j8yme3vszlrg.svg)

```mathematica
crossSectionTotal = 2*Pi*Integrate[diffXSection2, {t, -s, 0}]
```

![1w3j66hnfojwz](img/1w3j66hnfojwz.svg)

Notice that up to the overall factor color factor 3 and the quark electric charge squared this result is identical to the total cross-section for the muon production in electron-positron annihilation.

```mathematica
crossSectionTotalQED = 4*Pi*(SMP["alpha_fs"]^2/3/s)
```

![0tyybs0txs9cb](img/0tyybs0txs9cb.svg)

Taking the ratio of the two gives us the famous R-ration prediction of the parton mode, where the summation over the quark flavors in front of the charge squared is understood

```mathematica
crossSectionTotal/crossSectionTotalQED
```

![1828vxlbcwm8x](img/1828vxlbcwm8x.svg)

```mathematica
quarkCharges = {eq[u | c | t] -> 2/3, eq[d | s | b] -> 
         -3^(-1)}; 
```

Depending on the available center of mass energy, we may not be able to produce all the existing
quark flavors. Below 3 GeV (roughly twice the mass of the charm quark) we have only up, down and strange quarks and the R-ratio is given by

```mathematica
Sum[3*eq[i]^2, {i, {u, d, s}}] /. quarkCharges
```

![05lusigoq8qf9](img/05lusigoq8qf9.svg)

At higher energies but below 9 GeV (roughly twice the mass of the bottom quark) we also have the 
contribution from the charm quark

```mathematica
Sum[3*eq[i]^2, {i, {u, d, s, c}}] /. quarkCharges
```

![0sythvh55ot7v](img/0sythvh55ot7v.svg)

At even higher energies the bottom quark must also be taken into account

```mathematica
Sum[3*eq[i]^2, {i, {u, d, s, c, b}}] /. quarkCharges
```

![0z5hiplizq3uk](img/0z5hiplizq3uk.svg)

At some point we finally reach sufficiently high energies to produce the top quark

```mathematica
Sum[3*eq[i]^2, {i, {u, d, s, c, b, t}}] /. quarkCharges
```

![097v8u62j7s2n](img/097v8u62j7s2n.svg)

## Check the final results

```mathematica
knownResults = {(6*(t^2 + u^2)*SMP["e"]^4*SMP["e_Q"]^2)/s^2, 
       (4*Pi*SMP["alpha_fs"]^2*SMP["e_Q"]^2)/s}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0], 
     crossSectionTotal}, knownResults, 
   Text -> {"\tCompare to CalcHEP and to Field, Applications of \
   Perturbative QCD, Eq. 2.1.15:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![0yxbl7g9w2uz1](img/0yxbl7g9w2uz1.svg)

![1dadttaznnpd3](img/1dadttaznnpd3.svg)

![1gbpucdl5cil5](img/1gbpucdl5cil5.svg)