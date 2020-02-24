---
title: Top quark decaying into a quark and a W boson
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Qt -> Qb W, EW, total decay rate, tree"; 
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
```

Enable CKM mixing

```mathematica
$CKM = True; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {F[3, {3}]} -> {F[4, {3}], -V[3]}, InsertionLevel -> 
         {Particles}]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0nheuowy9ffrq](img/0nheuowy9ffrq.svg)

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

![1j6to5sun2aci](img/1j6to5sun2aci.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]
SP[p] = SMP["m_t"]^2; 
SP[k1] = SMP["m_b"]^2; 
SP[k2] = SMP["m_W"]^2; 
SP[k1, k2] = Simplify[(SP[p] - SP[k1] - SP[k2])/2]; 
SP[p, k1] = Simplify[ExpandScalarProduct[SP[k1 + k2, k1]]]; 
SP[p, k2] = Simplify[ExpandScalarProduct[SP[k1 + k2, k2]]]; 
```

## Square the amplitude

We average over the polarizations of the top quark, hence the additional factor 1/2

```mathematica
ampSquared[0] = Simplify[(DoPolarizationSums[#1, k2] & )[
       DiracSimplify[(FermionSpinSum[#1, ExtraFactor -> 1/2] & )[
           SUNSimplify[amp[0]*ComplexConjugate[amp[0]]]]]]]
```

![15d4p10n0d1lf](img/15d4p10n0d1lf.svg)

## Total decay rate

```mathematica
phaseSpacePrefactor[m1_, m2_, M_] := (1/(16*Pi*M))*
       Sqrt[1 - (m1 + m2)^2/M^2]*Sqrt[1 - (m1 - m2)^2/M^2]; 
```

```mathematica
totalDecayRate = 
   (#1 /. Sqrt[x_]*Sqrt[y_] :> Sqrt[ExpandAll[x*y]] & )[
     Simplify[phaseSpacePrefactor[SMP["m_b"], SMP["m_W"], 
           SMP["m_t"]]*ampSquared[0]]]
```

![1p5gw9esbobil](img/1p5gw9esbobil.svg)

## Check the final results

```mathematica
knownResults = {SMP["m_t"]^3*
         ((CA*SMP["G_F"]*Sqrt[((SMP["m_b"] - SMP["m_t"] - 
                         SMP["m_W"])*(SMP["m_b"] + SMP["m_t"] - 
              SMP["m_W"])*
                      (SMP["m_b"] - SMP["m_t"] + SMP["m_W"])*
                      (SMP["m_b"] + SMP["m_t"] + SMP["m_W"]))/
          SMP["m_t"]^4]*
               ((1 - SMP["m_b"]^2/SMP["m_t"]^2)^2 + 
                  (SMP["m_W"]^2/SMP["m_t"]^2)*(1 + SMP["m_b"]^2/
                         SMP["m_t"]^2) - 
          2*(SMP["m_W"]^4/SMP["m_t"]^4))*
               SMP["V_tb", -I]*SMP["V_tb", I])/(8*Sqrt[2]*Pi))}; 
FCCompareResults[{totalDecayRate}, knownResults, 
     Text -> {"\tCompare to Grozin, Using REDUCE in High Energy \
    Physics, Chapter 5.2:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 3], 0.001], 
     " s."]; 
```

![13887nuesqmlc](img/13887nuesqmlc.svg)

![04o7mzv7lvpty](img/04o7mzv7lvpty.svg)