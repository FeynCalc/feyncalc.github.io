---
title: Muon production from the decay of a virtual photon
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga^* -> Mu Amu, QED, total decay rate, tree"; 
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
MakeBoxes[p, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(p\), \(1\)]\)"; 
MakeBoxes[k1, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(1\)]\)"; 
MakeBoxes[k2, TraditionalForm] := 
     "\!\(\*SubscriptBox[\(k\), \(2\)]\)"; 
```

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[1]} -> {F[2, {2}], -F[2, {2}]}, InsertionLevel -> 
         {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {1, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {256, 256}]; 
```

![1i70msj6dzd1m](img/1i70msj6dzd1m.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True]
```

![0b2tulm64hhhb](img/0b2tulm64hhhb.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[k1] = SMP["m_mu"]^2; 
SP[k2] = SMP["m_mu"]^2; 
SP[k1, k2] = (QQ - SP[k1] - SP[k2])/2; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (DoPolarizationSums[#1, p, 0, VirtualBoson -> True] & )[
       DiracSimplify[FermionSpinSum[FeynAmpDenominatorExplicit[
             amp[0]*ComplexConjugate[amp[0]]]]]]]
```

![0qy3sqryncq7d](img/0qy3sqryncq7d.svg)

```mathematica
ampSquaredMassless[0] = (#1 /. {SMP["m_mu"] -> 0} & )[
     ampSquared[0]]
```

![1twril0tl5jql](img/1twril0tl5jql.svg)

## Total decay rate

The differential decay rate  d Gamma/ d Omega is given by

```mathematica
prefac = ExpandScalarProduct[(1/(64*Pi^2))*
       (1/Sqrt[SP[k1 + k2]])]
```

![00kp61gss62gf](img/00kp61gss62gf.svg)

```mathematica
diffDecayRate = prefac*ampSquaredMassless[0] /. 
     SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"]
```

![0p4h20kjccn6a](img/0p4h20kjccn6a.svg)

The total decay-rate

```mathematica
decayRateTotal = 4*Pi*diffDecayRate
```

![0ukqmx3prsi1k](img/0ukqmx3prsi1k.svg)

## Check the final results

```mathematica
knownResults = {SMP["alpha_fs"]*Sqrt[QQ]}; 
FCCompareResults[{decayRateTotal}, knownResults, 
     Text -> {"\tCompare to Field, Applications of Perturbative \
    QCD, Eq. 2.1.29", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![03dpogus4xfx5](img/03dpogus4xfx5.svg)

![1cn36yvrop07a](img/1cn36yvrop07a.svg)
