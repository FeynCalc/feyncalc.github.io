---
title: Quark-antiquark production from the decay of a virtual photon
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga^* -> Q Qbar, QCD, total decay rate, tree"; 
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

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 2], 
       {V[1]} -> {F[3, {1}], -F[3, {1}]}, InsertionLevel -> 
         {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![0un7xask8cval](img/0un7xask8cval.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     DropSumOver -> True, Prefactor -> (3/2)*SMP["e_Q"], 
     FinalSubstitutions -> {SMP["m_u"] -> SMP["m_q"]}]
```

![1h54j1juppzh6](img/1h54j1juppzh6.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[k1] = SMP["m_q"]^2; 
SP[k2] = SMP["m_q"]^2; 
SP[k1, k2] = (QQ - SP[k1] - SP[k2])/2; 
```

## Square the amplitude

```mathematica
ampSquared[0] = DiracSimplify[
     (DoPolarizationSums[#1, p, 0, VirtualBoson -> True] & )[
       FermionSpinSum[(SUNSimplify[#1, Explicit -> True, 
                SUNNToCACF -> False] & )[FeynAmpDenominatorExplicit[
             amp[0]*ComplexConjugate[amp[0]]]]]]]
```

![03w1tbu2rgk9c](img/03w1tbu2rgk9c.svg)

```mathematica
ampSquaredMassless[0] = (#1 /. {SMP["m_q"] -> 0} & )[
     ampSquared[0]]
```

![0vxwulqulag3r](img/0vxwulqulag3r.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![0phwv9fnm1ikq](img/0phwv9fnm1ikq.svg)

## Total decay rate

The differential decay rate  d Gamma/ d Omega is given by

```mathematica
prefac = ExpandScalarProduct[(1/(64*Pi^2))*
       (1/Sqrt[SP[k1 + k2]])]
```

![00kp61gss62gf](img/00kp61gss62gf.svg)

```mathematica
diffDecayRate = prefac*ampSquaredMasslessSUNN3[0] /. 
     SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"]
```

![1e6okob0w612p](img/1e6okob0w612p.svg)

The total decay-rate

```mathematica
decayRateTotal = 4*Pi*diffDecayRate
```

![1jmhjny5b0gb3](img/1jmhjny5b0gb3.svg)

Notice that up to the overall color factor 3 and the quark electric charge squared this result is identical to the total decay rate of a virtual photon into a muon-antimuon pair

```mathematica
decayRateTotalQED = SMP["alpha_fs"]*Sqrt[QQ]
```

![00fx3ldmu0kjp](img/00fx3ldmu0kjp.svg)

Taking the ration of the two gives us the famous R-ration prediction of the parton mode, where the summation over the quark flavors in front of the charge squared is understood

```mathematica
decayRateTotal/decayRateTotalQED
```

![00zzsx5jxuc34](img/00zzsx5jxuc34.svg)

## Check the final results

```mathematica
knownResults = {3*SMP["alpha_fs"]*SMP["e_Q"]^2*Sqrt[QQ]}; 
FCCompareResults[{decayRateTotal}, knownResults, 
     Text -> {"\tCompare to Field, Applications of Perturbative \
    QCD, Eq. 2.1.30", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0izn324hnkxa4](img/0izn324hnkxa4.svg)

![0xlm50lxix94p](img/0xlm50lxix94p.svg)