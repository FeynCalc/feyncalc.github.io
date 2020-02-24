---
title: Adler-Bell-Jackiw anomaly in QED
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "Pi -> Ga Ga, QED, axial current, 1-loop"; 
If[$FrontEnd === Null, $FeynCalcStartupMessages = False; 
      Print[description]; ]; 
If[$Notebooks === False, $FeynCalcStartupMessages = False]; 
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

## Obtain the amplitude

Nicer typesetting

```mathematica
MakeBoxes[mu, TraditionalForm] := "\[Mu]"; 
MakeBoxes[nu, TraditionalForm] := "\[Nu]"; 
MakeBoxes[la, TraditionalForm] := "\[Lambda]"; 
```

According to Peskin and Schroeder (Ch 19.2), the amplitude for the first triangle diagram reads

```mathematica
amp1[0] = Explicit[-(((-I)*SMP["e"])^2*
          DiracTrace[GAD[mu] . GA[5] . QuarkPropagator[l - k] . 
              GAD[la] . QuarkPropagator[l] . GAD[nu] . 
              QuarkPropagator[l + p]])]
```

![0x1riy0vfghrt](img/0x1riy0vfghrt.svg)

And the second one follows from the first by interchanging k with p and la with nu

```mathematica
amp2[0] = amp1[0] /. {k -> p, p -> k, la -> nu, nu -> la}
```

![1ipvi3323d4nx](img/1ipvi3323d4nx.svg)

## Calculate the amplitude

Contracting both amplitudes with I*(k+p)^mu we can check the non-conservation of the axial current.

```mathematica
amp[0] = Contract[I*FVD[k + p, mu]*(amp1[0] + amp2[0])]
```

![0qkpnfkuy8gu8](img/0qkpnfkuy8gu8.svg)

For this calculation it is crucial to use a correct scheme for gamma^5. As in the book, we use the 
Breitenlohner-Maison-t'Hooft-Veltman prescription.

```mathematica
FCSetDiracGammaScheme["BMHV"]; 
amp[1] = TID[amp[0], l, ToPaVe -> True]
```

![06yj54kyu4l6i](img/06yj54kyu4l6i.svg)

```mathematica
FCClearScalarProducts[]; 
Momentum[k, D | D - 4] = Momentum[k]; 
Momentum[p, D | D - 4] = Momentum[p]; 
```

The explicit values for the PaVe functions B0 and C0 can be obtained e.g. from H. Patel's Package-X. 
Here we just insert the known results. The C0 function is finite here, so because of the prefactor (D-4) it 
gives no contribution in the D->4 limit.

```mathematica
amp[2] = Collect2[amp[1], {B0, C0}] //. 
     {B0[FCI[SP[p_, p_]], 0, 0] :> 1/(16*Epsilon*Pi^4) - 
           (-2 + EulerGamma)/(16*Pi^4) + 
           Log[-((4*Pi*ScaleMu^2)/Pair[Momentum[p], Momentum[p]])]/
             (16*Pi^4), B0[FCI[SP[p, p] + 2*SP[p, k] + SP[k, k]], 0, 
           0] :> B0[FCI[SP[k + p, k + p]], 0, 0], 
       (D - 4)*ExpandScalarProduct[C0[SP[k], SP[p], SP[k + p], 0, 
               0, 0]] -> 0}
```

![0xj4q7zdlm5jp](img/0xj4q7zdlm5jp.svg)

Now we insert the explicit values, convert the external momenta to 4 dimensions and expand in Epsilon

```mathematica
amp[3] = Normal[(Series[#1, {Epsilon, 0, 0}] & )[
       (FCReplaceD[#1, D -> 4 - 2*Epsilon] & )[amp[2]]]]
```

![0dl636wc83mau](img/0dl636wc83mau.svg)

The result should be twice Eq. 19.59 in Peskin and Schroeder

## Check the final results

```mathematica
knownResult = Contract[2*((SMP["e"]^2/(4*Pi^2))*
            LC[al, la, be, nu]*FV[k, al]*FV[p, be])]; 
FCCompareResults[amp[3], knownResult, 
     Text -> {"\tCompare to Peskin and Schroeder, An Introduction \
    to QFT, Eq 19.59:", "CORRECT.", "WRONG!"}, 
     Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0kdw4zh9irj2d](img/0kdw4zh9irj2d.svg)

![1etajoox1jga1](img/1etajoox1jga1.svg)