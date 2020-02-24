---
title: Electron's g-2 in QED
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "El -> Ga El, QED, F2(0) form factor, 1-loop"; 
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
MakeBoxes[mu, TraditionalForm] := "\[Mu]"; 
MakeBoxes[p1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(1\)]\)"; 
MakeBoxes[p2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[1, 1 -> 2, 
         ExcludeTopologies -> {Tadpoles, WFCorrections}], 
       {F[2, {1}]} -> {V[1], F[2, {1}]}, InsertionLevel -> 
         {Particles}, ExcludeParticles -> {S[_], V[2 | 3], 
           (S | U)[_], F[3 | 4], F[2, {2 | 3}]}]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![04h1g9ef17ysa](img/04h1g9ef17ysa.svg)

## Obtain the amplitude

The 1/(2Pi)^D prefactor is implicit. We need to replace e with -e to be compatible
with the convnention D^mu = d^mu + ie A^mu

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags, PreFactor -> 1], 
         IncomingMomenta -> {p1}, OutgoingMomenta -> {k, p2}, 
         LorentzIndexNames -> {mu}, LoopMomenta -> {q}, 
         UndoChiralSplittings -> True, ChangeDimension -> D, 
         List -> False, SMP -> True, FinalSubstitutions -> 
           {SMP["e"] -> -SMP["e"]}] /. k -> p1 - p2 /. q -> q + p1
```

![1xtmgipdb4lha](img/1xtmgipdb4lha.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
ME = SMP["m_e"]; 
ScalarProduct[p1, p1] = ME^2; 
ScalarProduct[p2, p2] = ME^2; 
ScalarProduct[k, k] = 0; 
ScalarProduct[p1, p2] = ME^2; 
```

## Calculate the amplitude

Amputate the polarization vector.

```mathematica
amp[1] = (#1 /. SMP["e"]^3 -> 4*Pi*SMP["e"]*SMP["alpha_fs"] & )[
     Contract[(#1 /. Pair[Momentum[Polarization[___], ___], 
                  ___] :> 1 & )[amp[0]]]]
```

![1veymtkz0l33c](img/1veymtkz0l33c.svg)

```mathematica
amp[2] = (Collect2[#1, Spinor] & )[DiracSimplify[
       TID[amp[1], q, ToPaVe -> True]]]
```

![0mk2tmho9r1oi](img/0mk2tmho9r1oi.svg)

To extract F2 (0) we need to look only at the piece proportional to (p1+p2)^mu. Thus we can drop the g^mu -piece

```mathematica
amp[3] = DotSimplify[(#1 /. FCI[GAD[mu]] :> 0 & )[amp[2]]]
```

![00ctjm5m8o0ex](img/00ctjm5m8o0ex.svg)

The explicit values for the PaVe functions C1, C11 and C12 can be obtained e.g. from H. Patel's Package-X. Here we just insert the known results.

```mathematica
amp[4] = amp[3] /. {PaVe[1, {ME^2, 0, ME^2}, {0, ME^2, ME^2}, 
           OptionsPattern[]] -> 1/(32*Pi^4*ME^2), 
       PaVe[1, 1, {ME^2, 0, ME^2}, {0, ME^2, ME^2}, 
           OptionsPattern[]] -> -(1/(96*Pi^4*ME^2)), 
       PaVe[1, 2, {ME^2, 0, ME^2}, {0, ME^2, ME^2}, 
           OptionsPattern[]] -> -(1/(192*Pi^4*ME^2))}
```

![0gxf6goytnowq](img/0gxf6goytnowq.svg)

As expected, F2 (0) is free of any divergences. So we can safely do the limit D ->4

```mathematica
amp[5] = (#1 /. D -> 4 & )[(ChangeDimension[#1, 4] & )[amp[4]]]
```

![0drd89jvwcbjp](img/0drd89jvwcbjp.svg)

We obtained $\frac{i e}{2 m_e} (p_1+p_2)^\mu F_2 (0) \bar{u}(p_2) u(p_1)$.
Dividing by the numerical prefactor and substituting $e^2 = 4\pi^2 \alpha$ yields F2(0)

```mathematica
f2[0] = (#1 /. {Spinor[__] . Spinor[__] :> 1, 
            FCI[FV[p1, _] + FV[p2, _]] :> 1} & )[
     amp[5]/((I*SMP["e"])/(2*ME))]
```

![0x5yttlk4rexg](img/0x5yttlk4rexg.svg)

## Check the final results

```mathematica
knownResult = AlphaFS/(2*Pi); 
FCCompareResults[f2[0], knownResult, 
     Text -> 
       {
         "\tCompare to J. Schwinger, Phys. Rev. 73, 416-417, 1948:"\
     , "CORRECT.", "WRONG!"}, Interrupt -> {Hold[Quit[1]], 
         Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![11nzba0k844uk](img/11nzba0k844uk.svg)

![12kzeb4ts42yh](img/12kzeb4ts42yh.svg)