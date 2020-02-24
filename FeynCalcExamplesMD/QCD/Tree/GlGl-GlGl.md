---
title: Gluon-gluon to gluon-gluon scattering
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "Gl Gl -> Gl Gl, QCD, matrix element squared, tree"; 
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
MakeBoxes[k3, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(3\)]\)"; 
MakeBoxes[k4, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(4\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 2 -> 2], 
       {V[5], V[5]} -> {V[5], V[5]}, InsertionLevel -> {Classes}, 
       Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1lgpgwezzcwv0](img/1lgpgwezzcwv0.svg)

![16s88hoxkfxol](img/16s88hoxkfxol.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {k1, k2}, OutgoingMomenta -> {k3, k4}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1, k2, k3, k4}, 
     List -> True, SMP -> True, Contract -> True, 
     DropSumOver -> True]
```

![1ezzubjstacjd](img/1ezzubjstacjd.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, k1, k2, -k3, -k4, 0, 0, 0, 0]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = (SUNSimplify[#1, Explicit -> True, 
            SUNNToCACF -> False] & )[FeynAmpDenominatorExplicit[
         (1/(SUNN^2 - 1)^2)*(amp[1]*ComplexConjugate[amp[0]])]]; 
```

```mathematica
polsums[x_, vec_, aux_, spinfac_] := 
   (FixedPoint[ReleaseHold, #1] & )[
     (DoPolarizationSums[#1, vec, aux, ExtraFactor -> 
              spinfac] & )[(Isolate[#1, {Polarization[vec, __]}] & )[
         (Collect2[#1, Pair[_, Momentum[Polarization[vec, 
                      __]]]] & )[x]]]]
```

```mathematica
ClearAll[re]; 
Table[
     Print["    calculating color factors in products of the \
   amplitudes ", i, " and ", j, " (CC), time = ", 
        Timing[re[i, j] = (SUNSimplify[#1, Explicit -> True, 
                     SUNNToCACF -> False] & )[
       FeynAmpDenominatorExplicit[
                  
        amp[0][[i]]*ComplexConjugate[amp[0]][[j]]]]][[1]]]; 
      re[i, j], {i, 4}, {j, i}]; 
```

![0u9nrcky1geul](img/0u9nrcky1geul.svg)

![1v15nyysckna8](img/1v15nyysckna8.svg)

![1d7dn9fn1bz8d](img/1d7dn9fn1bz8d.svg)

![10702l4ci8wi9](img/10702l4ci8wi9.svg)

![0zegbfpxyaill](img/0zegbfpxyaill.svg)

![19mg85qheaeba](img/19mg85qheaeba.svg)

![0c66mb9utvoid](img/0c66mb9utvoid.svg)

![01igr5axax858](img/01igr5axax858.svg)

![0ei6evk4tm8vv](img/0ei6evk4tm8vv.svg)

![1443azlempyb9](img/1443azlempyb9.svg)

```mathematica
ClearAll[pre]; 
Table[Print["    calculating product of the amplitudes ", i, 
        " and ", j, " (CC), time = ", 
        Timing[pre[i, j] = Simplify[(polsums[#1, k4, k3, 1] & )[
                  (polsums[#1, k3, k4, 1] & )[
                    (polsums[#1, k2, k1, 1/2] & )[
                      (polsums[#1, k1, k2, 1/2] & )[
           re[i, j]]]]]]][[1]]]; 
      pre[i, j], {i, 4}, {j, i}]; 
```

![0c3bbz646r8es](img/0c3bbz646r8es.svg)

![0ihn8lhyixbyr](img/0ihn8lhyixbyr.svg)

![1qvdxufol3qm7](img/1qvdxufol3qm7.svg)

![03wwfptmebht2](img/03wwfptmebht2.svg)

![1nwn72j4zkp6b](img/1nwn72j4zkp6b.svg)

![099sa9jtix2bz](img/099sa9jtix2bz.svg)

![0hil6zeqcrsd5](img/0hil6zeqcrsd5.svg)

![1u3kcuxtcrefz](img/1u3kcuxtcrefz.svg)

![18aslnmjr427k](img/18aslnmjr427k.svg)

![15tmz3dygeo0p](img/15tmz3dygeo0p.svg)

```mathematica
fpre[i_, j_] := pre[i, j] /; i >= j; 
fpre[i_, j_] := ComplexConjugate[pre[j, i]] /; i < j; 
ampSquared[0] = Simplify[(TrickMandelstam[#1, {s, t, u, 0}] & )[
       (1/(SUNN^2 - 1)^2)*Sum[fpre[i, j], {i, 1, 4}, {j, 1, 4}]]]
```

![0shhdot80hkw2](img/0shhdot80hkw2.svg)

```mathematica
ampSquaredSUNN3[0] = ampSquared[0] /. SUNN -> 3
```

![0nn780gxmrf37](img/0nn780gxmrf37.svg)

```mathematica
ampSquaredMassless[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     (#1 /. {SMP["m_u"] -> 0} & )[ampSquared[0]]]
```

![08u82reaxtbv9](img/08u82reaxtbv9.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![17lk0hfhy9p53](img/17lk0hfhy9p53.svg)

## Check the final results

```mathematica
knownResults = {(9/2)*SMP["g_s"]^4*(3 - t*(u/s^2) - s*(u/t^2) - 
            s*(t/u^2))}; 
FCCompareResults[{ampSquaredMasslessSUNN3[0]}, {knownResults}, 
   Text -> {"\tCompare to Ellis, Stirling and Weber, QCD and \
   Collider Physics, Table 7.1:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}, 
   Factoring -> Function[x, Simplify[TrickMandelstam[x, 
           {s, t, u, 0}]]]]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![19r5z3aom1xl6](img/19r5z3aom1xl6.svg)

![06nhs1ngsr7ts](img/06nhs1ngsr7ts.svg)

![0crmy8e5ou1g3](img/0crmy8e5ou1g3.svg)