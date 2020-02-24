---
title: W decaying into a quark and an antiquark of different flavors
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "W -> Qi Qjbar, EW, total decay rate, tree"; 
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

Enable CKM mixing

```mathematica
$CKM = True; 
```

Nicer typesetting

```mathematica
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[3]} -> {-F[3, {1}], F[4, {1}]}, InsertionLevel -> 
         {Particles}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![07lhlnsgwzovw](img/07lhlnsgwzovw.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     ChangeDimension -> 4, List -> False, SMP -> True, 
     Contract -> True, DropSumOver -> True, 
     FinalSubstitutions -> 
       {SMP["e"] -> Sqrt[(8/Sqrt[2])*SMP["G_F"]*SMP["m_W"]^2*
               SMP["sin_W"]^2]}]
```

![0iq8nee7efz1a](img/0iq8nee7efz1a.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]
SP[p] = SMP["m_W"]^2; 
SP[k1] = SMP["m_u"]^2; 
SP[k2] = SMP["m_d"]^2; 
SP[k1, k2] = (SMP["m_W"]^2 - SMP["m_u"]^2 - SMP["m_d"]^2)/2; 
SP[p, k1] = Simplify[ExpandScalarProduct[SP[k1 + k2, k1]]]; 
SP[p, k2] = Simplify[ExpandScalarProduct[SP[k1 + k2, k2]]]; 
```

## Square the amplitude

We average over the polarizations of the W-boson, hence the additional factor 1/3.

```mathematica
ampSquared[0] = Simplify[
     (DoPolarizationSums[#1, p, ExtraFactor -> 1/3] & )[
       DiracSimplify[FermionSpinSum[SUNSimplify[
             amp[0]*ComplexConjugate[amp[0]]]]]]]
```

![0bdvgj9j6g66s](img/0bdvgj9j6g66s.svg)

## Total decay rate

```mathematica
phaseSpacePrefactor[m1_, m2_, M_] := (1/(16*Pi*M))*
       Sqrt[1 - (m1 + m2)^2/M^2]*Sqrt[1 - (m1 - m2)^2/M^2]; 
```

```mathematica
totalDecayRate = Factor2[
     (#1 /. Sqrt[x_]*Sqrt[y_] :> Sqrt[ExpandAll[x*y]] & )[
       Simplify[phaseSpacePrefactor[SMP["m_u"], SMP["m_d"], 
             SMP["m_W"]]*ampSquared[0]]]]
```

![0su9wf1qcq1lf](img/0su9wf1qcq1lf.svg)

## Check the final results

```mathematica
knownResults = {SMP["m_W"]^3*
         ((CA*SMP["G_F"]*Sqrt[((SMP["m_d"] - SMP["m_u"] - 
                         SMP["m_W"])*(SMP["m_d"] + SMP["m_u"] - 
              SMP["m_W"])*
                      (SMP["m_d"] - SMP["m_u"] + SMP["m_W"])*
                      (SMP["m_d"] + SMP["m_u"] + SMP["m_W"]))/
          SMP["m_W"]^4]*
               (1 - (SMP["m_u"]^2 + SMP["m_d"]^2)/(2*SMP["m_W"]^2) - 
                  (SMP["m_u"]^2 - SMP["m_d"]^2)^2/(2*SMP["m_W"]^4))*
               SMP["V_ud", -I]*SMP["V_ud", I])/(6*Sqrt[2]*Pi))}; 
FCCompareResults[{totalDecayRate}, knownResults, 
   Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
   Physics, Chapter 5.2:", "CORRECT.", "WRONG!"}, 
   Interrupt -> {Hold[Quit[1]], Automatic}]
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![13887nuesqmlc](img/13887nuesqmlc.svg)

![00lgxwvusjos5](img/00lgxwvusjos5.svg)

![0v9ssqqkf738z](img/0v9ssqqkf738z.svg)