---
title: Virtual photon-gluon scattering to a quark-antiquark pair
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Ga Gl -> Q Qbar, QCD, matrix element squared, tree"; 
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
       {V[1], V[5]} -> {F[3, {1}], -F[3, {1}]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0escq6x6d0p70](img/0escq6x6d0p70.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {p2}, List -> False, 
     SMP -> True, Contract -> True, DropSumOver -> True, 
     Prefactor -> (3/2)*SMP["e_Q"]]
```

![1qvjhly7zf1y1](img/1qvjhly7zf1y1.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, qQ, 0, SMP["m_u"], 
     SMP["m_u"]]; 
```

## Square the amplitude

Now come the usual steps, but with some special features. We do not average over the polarizations of the virtual photon, but use the gauge trick for the sum over its polarizations. In this case the sum goes over all 4 unphysical polarizations,  not just 2.

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_u"]^2 + qQ^2}] & )[
       (DoPolarizationSums[#1, p2, k1] & )[
         (DoPolarizationSums[#1, p1, 0, VirtualBoson -> True, 
                GaugeTrickN -> 4] & )[DiracSimplify[
             (FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
               (SUNSimplify[#1, Explicit -> True, SUNNToCACF -> 
                        False] & )[FeynAmpDenominatorExplicit[
                   (1/(SUNN^2 - 1))*(amp[0]*ComplexConjugate[
                          amp[0]])]]]]]]]]
```

![00vo3yh5tlxcu](img/00vo3yh5tlxcu.svg)

```mathematica
ampSquaredMassless[0] = 
   (TrickMandelstam[#1, {s, t, u, qQ^2}] & )[
     (#1 /. {SMP["m_u"] -> 0} & )[ampSquared[0]]]
```

![0avsj8vnlqh3l](img/0avsj8vnlqh3l.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = 
   Simplify[ampSquaredMassless[0] /. SUNN -> 3 /. qQ -> I*Q]
```

![1vybh7f4bwszw](img/1vybh7f4bwszw.svg)

## Check the final results

```mathematica
knownResults = {SMP["e"]^2*SMP["e_Q"]^2*SMP["g_s"]^2*2*
         (u/t + t/u + 2*Q^2*((u + t + Q^2)/(t*u)))}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0]}, {knownResults}, 
     Text -> {"\tCheck with R. Field, Applications of \
    Perturbative QCD, Eq 4.3.20:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![06vfk3l7r16os](img/06vfk3l7r16os.svg)

![0z4v58vogntzj](img/0z4v58vogntzj.svg)