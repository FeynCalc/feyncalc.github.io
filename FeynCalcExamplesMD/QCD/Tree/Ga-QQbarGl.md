---
title: Quark-antiquark production with a real gluon emission from the decay of a virtual photon
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Ga^* -> Q Qbar Gl, QCD, total decay rate, tree"; 
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
```

```mathematica
diags = InsertFields[CreateTopologies[0, 1 -> 3], 
       {V[1]} -> {F[3, {1}], -F[3, {1}], V[5]}, 
       InsertionLevel -> {Classes}, Model -> "SMQCD"]; 
Paint[diags, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {512, 256}]; 
```

![103mav02fxn0e](img/103mav02fxn0e.svg)

## Obtain the amplitude

```mathematica
amp[0] = FCFAConvert[CreateFeynAmp[diags], 
     IncomingMomenta -> {p}, OutgoingMomenta -> {k1, k2, k3}, 
     UndoChiralSplittings -> True, ChangeDimension -> 4, 
     List -> False, SMP -> True, Contract -> True, 
     DropSumOver -> True, Prefactor -> (3/2)*SMP["e_Q"], 
     FinalSubstitutions -> {SMP["m_u"] -> SMP["m_q"]}]
```

![0lxrfgo8qfhvo](img/0lxrfgo8qfhvo.svg)

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SP[k1] = SMP["m_q"]^2; 
SP[k2] = SMP["m_q"]^2; 
SP[k3] = 0; 
SP[k1, k2] = (QQ/2)*(1 - x3); 
SP[k1, k3] = (QQ/2)*(1 - x2); 
SP[k2, k3] = (QQ/2)*(1 - x1); 
```

## Square the amplitude

```mathematica
ampSquared[0] = Simplify[FeynAmpDenominatorExplicit[
       DiracSimplify[FermionSpinSum[
           (DoPolarizationSums[#1, k3, 0, VirtualBoson -> True] & )[
             (DoPolarizationSums[#1, p, 0, VirtualBoson -> True] & )[
               SUNSimplify[amp[0]*ComplexConjugate[amp[0]]]]]]]]]
```

![0lod276p66z3o](img/0lod276p66z3o.svg)

```mathematica
ampSquaredMassless[0] = 
   (SUNSimplify[#1, SUNNToCACF -> False] & )[
     Simplify[(#1 /. {SMP["m_q"] -> 0, x3 -> 2 - x1 - x2, 
                SMP["e"]^2 -> 4*Pi*SMP["alpha_fs"], SMP["g_s"]^2 -> 
                  4*Pi*SMP["alpha_s"]} & )[ampSquared[0]]]]
```

![0b3ayqko3gfun](img/0b3ayqko3gfun.svg)

```mathematica
ampSquaredMasslessSUNN3[0] = ampSquaredMassless[0] /. SUNN -> 3
```

![0pjoo8cnynex8](img/0pjoo8cnynex8.svg)

## Total decay rate

```mathematica
pref = (QQ/(128*Pi^3))*(1/(2*Sqrt[QQ]))
```

![0tbtcsxcxkul4](img/0tbtcsxcxkul4.svg)

```mathematica
normBorn = 3*SMP["alpha_fs"]*SMP["e_Q"]^2*Sqrt[QQ]
```

![0dmeka2j5kshk](img/0dmeka2j5kshk.svg)

Differential cross-section normalized w.r.t to the Born cross-section 1/sigma_0 d sigma / (d x1 d x2)

```mathematica
normDiffCrossSection = ampSquaredMasslessSUNN3[0]*
     (pref/normBorn)
```

![07e86z3j5j85s](img/07e86z3j5j85s.svg)

This integral is divergent for x1->1 and x2->1. The source of these divergences are infrared (when the gluon energy
approaches 0)  and collinear (when the gluon and quark become collinear) singularities.

```mathematica
If[$FrontEnd =!= Null, Plot3D[normDiffCrossSection /. 
       SMP["alpha_s"] -> 1, {x1, 0, 1}, {x2, 0, 1}]]
```

![058o78m3r81pt](img/058o78m3r81pt.svg)

Introducing a regulator beta=m^2/Q^2 to enforce that the Mandelstam variables s and t are always larger than m^2 gives

```mathematica
normDiffCrossSection
```

![0zoe546eflk5g](img/0zoe546eflk5g.svg)

```mathematica
tmpIntegral = Integrate[normDiffCrossSection, 
     {x2, 1 - x1, 1 - beta}, Assumptions -> 
       {beta < x1, beta > 0, x1 >= 0, x1 <= 1}]
```

![0gcm8ilacol0k](img/0gcm8ilacol0k.svg)

```mathematica
integralReg = ConditionalExpression[
     ((5 - 10*beta - 4*(3 + (-4 + beta)*beta + (2*I)*Pi)*
               ArcTanh[1 - 2*beta] + 2*Log[1 - beta]*
               Log[(1 - beta)/beta^2] + 2*Log[beta]^2 + 
             4*PolyLog[2, (1 - beta)^(-1)] - 
       4*PolyLog[2, beta^(-1)])*
          SMP["alpha_s"])/(3*Pi), beta < 1/2]
```

![07g6o4t98cxng](img/07g6o4t98cxng.svg)

Expanding around beta=0 we obtain

```mathematica
integralRegExpanded = 
   Normal[Series[Simplify[Normal[integralReg]], {beta, 0, 0}, 
       Assumptions -> beta > 0]]
```

![1nau7mdrswjrm](img/1nau7mdrswjrm.svg)

Factoring out the Born cross-section we arrive to

```mathematica
integralRegExpandedFinal = Collect2[integralRegExpanded, Log, 
     Pi, FCFactorOut -> (2/(3*Pi))*SMP["alpha_s"]]
```

![0qvg0vaiialbw](img/0qvg0vaiialbw.svg)

To get rid of the singularities we must also include the virtual contributions to the cross-section!

## Check the final results

```mathematica
knownResults = {(2*(x1^2 + x2^2)*SMP["alpha_s"])/
         (3*Pi*(-1 + x1)*(-1 + x2))}; 
FCCompareResults[{normDiffCrossSection}, knownResults, 
     Text -> {"\tCompare to Field, Applications of Perturbative \
    QCD, Eq. 2.3.32", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![18jx38i5o1s54](img/18jx38i5o1s54.svg)

![08c4ym26lky5i](img/08c4ym26lky5i.svg)