---
title: Pair annihilation
---


```mathematica
Get["M2MD`"]
```

```mathematica
$toPlainText[cell_, ___] := StringJoin["---\ntitle: ", 
       M2MD`Private`BoxesToString[cell, "PlainText"], "\n---\n"]; 
```

```mathematica
MDExport["/home/vs/Downloads/outputMD/test.md", 
   Cell["Pair annihilation", "Title"], "CellStyleRules" -> 
     Association["Title" -> {"Text", $toPlainText}]]
```

![0p8zp61qbvcnd](img/0p8zp61qbvcnd.svg)

## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = 
     "El Ael -> Ga Ga, QED, matrix element squared, tree"; 
If[$FrontEnd === Null, $FeynCalcStartupMessages = False; 
      Print[description]; ]; 
If[$Notebooks === False, $FeynCalcStartupMessages = False]; 
$LoadAddOns = {"FeynArts"}; 
Get["FeynCalc`"]
$FAVerbose = 0; 
FCCheckVersion[9, 3, 0]; 
```

![0v0ipz351wwvy](img/0v0ipz351wwvy.svg)

![05ccutuzh8htk](img/05ccutuzh8htk.svg)

![1qzul9wezozfq](img/1qzul9wezozfq.svg)

![1cxfjk19kq9cn](img/1cxfjk19kq9cn.svg)

![1o49s8tpy7lgl](img/1o49s8tpy7lgl.svg)

![0vwaa9iez5qno](img/0vwaa9iez5qno.svg)

![14v35110v9xrm](img/14v35110v9xrm.svg)

![07kerhv4s5gyz](img/07kerhv4s5gyz.svg)

![1p7hfcd6vyiug](img/1p7hfcd6vyiug.svg)

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
       {F[2, {1}], -F[2, {1}]} -> {V[1], V[1]}, 
       InsertionLevel -> {Classes}, Restrictions -> QEDOnly]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![1p4as3biuav8p](img/1p4as3biuav8p.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     TransversePolarizationVectors -> {k1, k2}, List -> False, 
     SMP -> True, Contract -> True]
```

![1nq28735quuvz](img/1nq28735quuvz.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"], 
     SMP["m_e"], 0, 0]; 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[
     (TrickMandelstam[#1, {s, t, u, 2*SMP["m_e"]^2}] & )[
       DiracSimplify[(FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           (DoPolarizationSums[#1, k2, 0] & )[
             (DoPolarizationSums[#1, k1, 0] & )[
               FeynAmpDenominatorExplicit[amp[0]*ComplexConjugate[
                     amp[0]]]]]]]]]
```

![0h3qjd444sjuk](img/0h3qjd444sjuk.svg)

## Check the final results

```mathematica
knownResult = 2*SMP["e"]^4*(SP[p1, k2]/SP[p1, k1] + 
          SP[p1, k1]/SP[p1, k2] + 2*SMP["m_e"]^2*
            (1/SP[p1, k1] + 1/SP[p1, k2]) - SMP["m_e"]^4*
            (1/SP[p1, k1] + 1/SP[p1, k2])^2); 
FCCompareResults[ampSquared[0], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eq 5.105:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![1eb7ypb3gcjnu](img/1eb7ypb3gcjnu.svg)

![17mbny8mgqo93](img/17mbny8mgqo93.svg)
