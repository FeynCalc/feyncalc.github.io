---
title: Quark-virtual photon scattering to a gluon and a quark
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Q Ga -> Gl Q, QCD, matrix element squared, tree"; 
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
       {F[3, {1}], V[1]} -> {V[5], F[3, {1}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0vezdere5pd0l](img/0vezdere5pd0l.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1}, List -> False, 
     SMP -> True, Contract -> True, DropSumOver -> True, 
     Prefactor -> (3/2)*SMP["e_Q"]]
```

![04zjajjvtstwm](img/04zjajjvtstwm.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_u"], qQ, 0, 
     SMP["m_u"]]; 
```

## Square the amplitude

Now come the usual steps, but with some special features. We do not average over the polarizations of the virtual photon, but use the gauge trick for the sum over its polarizations. In this case the sum goes over all 4 unphysical polarizations,  not just 2.

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2 + qQ^2}] & )[
       (DoPolarizationSums[#1, k1, p2] & )[
         (DoPolarizationSums[#1, p2, 0, VirtualBoson -> True, 
                GaugeTrickN -> 4] & )[DiracSimplify[
             (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
               (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                        False] & )[FeynAmpDenominatorExplicit[
                   (1/SUNN)*(amp[0]*ComplexConjugate[amp[0]])]]]]]]]]
```

![06s4kxyc1uy06](img/06s4kxyc1uy06.svg)

```mathematica
ampSquaredMassless[0] = 
   (TrickMandelstam[#1, {s, t, u, qQ^2}] & )[
     (#1 /. {SMP["m_u"] -> 0} & )[ampSquared[0]]]
```

![1odz8qx2mpxmk](img/1odz8qx2mpxmk.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = 
   Simplify[ampSquaredMassless[0] /. SUNN -> 3 /. 
         u -> qQ^2 - s - t /. qQ -> I*Q]
```

![0w0u08owcfjxp](img/0w0u08owcfjxp.svg)

## Check the final results

```mathematica
knownResults = {(16/3)*SMP["e"]^2*SMP["e_Q"]^2*SMP["g_s"]^2*
         (-t/s - s/t - 2*Q^2*((s + t + Q^2)/(t*s)))}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0]}, {knownResults}, 
     Text -> {"\tCheck with R. Field, Applications of \
    Perturbative QCD, Eq 4.3.10:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![0xpieivuf0b40](img/0xpieivuf0b40.svg)

![0pjbrwiewj99o](img/0pjbrwiewj99o.svg)