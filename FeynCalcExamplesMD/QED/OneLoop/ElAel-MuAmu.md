---
title: Muon production
---


## Load FeynCalc and the necessary add-ons or other packages

```mathematica
description = "El Ael -> Mu Amu, QED, Born-virtual, 1-loop"; 
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
diagsTree = InsertFields[CreateTopologies[0, 2 -> 2, 
         ExcludeTopologies -> {Tadpoles, WFCorrections}], 
       {F[2, {1}], -F[2, {1}]} -> {F[2, {2}], -F[2, {2}]}, 
       InsertionLevel -> {Particles}, Restrictions -> QEDOnly, 
       ExcludeParticles -> {F[1 | 3 | 4, _], F[2, {3}]}]; 
Paint[diagsTree, ColumnsXRows -> {2, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {1024, 256}]; 
```

![0i97a6w6ho48d](img/0i97a6w6ho48d.svg)

```mathematica
diagsLoop = InsertFields[CreateTopologies[1, 2 -> 2, 
         ExcludeTopologies -> {Tadpoles, WFCorrections}], 
       {F[2, {1}], -F[2, {1}]} -> {F[2, {2}], -F[2, {2}]}, 
       InsertionLevel -> {Particles}, Restrictions -> QEDOnly, 
       ExcludeParticles -> {F[1 | 3 | 4, _], F[2, {3}]}]; 
Paint[DiagramExtract[diagsLoop, ((1) ..)*5], 
     ColumnsXRows -> {5, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {1024, 196}]; 
```

![0839xrz4y8nl5](img/0839xrz4y8nl5.svg)

```mathematica
diagsLoopCT = InsertFields[CreateCTTopologies[1, 2 -> 2, 
         ExcludeTopologies -> {Tadpoles, WFCorrectionCTs}], 
       {F[2, {1}], -F[2, {1}]} -> {F[2, {2}], -F[2, {2}]}, 
       InsertionLevel -> {Particles}, Restrictions -> QEDOnly, 
       ExcludeParticles -> {F[1 | 3 | 4, _], F[2, {3}]}]; 
Paint[diagsLoopCT, ColumnsXRows -> {3, 1}, Numbering -> Simple, 
     SheetHeader -> None, ImageSize -> {1024, 196}]; 
```

![0f6hb1pw9zsm8](img/0f6hb1pw9zsm8.svg)

## Obtain the amplitudes

```mathematica
ampLoopCT[0] = FCFAConvert[CreateFeynAmp[diagsLoopCT, 
           Truncated -> False, PreFactor -> 1] //. 
         {(h : dZfL1 | dZfR1)[z__] :> dZf1[z], 
           Conjugate[(h : dZfL1 | dZfR1)[z__]] :> dZf1[z], 
           dZZA1 -> 0}, IncomingMomenta -> {p1, p2}, 
       OutgoingMomenta -> {k1, k2}, LoopMomenta -> {l}, 
       ChangeDimension -> D, DropSumOver -> True, 
       UndoChiralSplittings -> True, SMP -> True, 
       FinalSubstitutions -> {SMP["m_e"] -> 0, SMP["m_mu"] -> 0}]; 
```

```mathematica
ampLoop[0] = FCFAConvert[CreateFeynAmp[DiagramExtract[
           diagsLoop, ((1) ..)*5], Truncated -> False, 
         PreFactor -> 1], IncomingMomenta -> {p1, p2}, 
       OutgoingMomenta -> {k1, k2}, LoopMomenta -> {q}, 
       ChangeDimension -> D, DropSumOver -> True, 
       UndoChiralSplittings -> True, SMP -> True, 
       FinalSubstitutions -> {SMP["m_e"] -> 0, SMP["m_mu"] -> 0}]; 
```

```mathematica
ampTree[0] = FCFAConvert[CreateFeynAmp[diagsTree, 
         Truncated -> False, PreFactor -> 1], 
       IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {k1, k2}, 
       ChangeDimension -> D, DropSumOver -> True, 
       UndoChiralSplittings -> True, SMP -> True, 
       FinalSubstitutions -> {SMP["m_e"] -> 0, SMP["m_mu"] -> 0}]; 
```

## Fix the kinematics

```mathematica
FCClearScalarProducts[]; 
SetMandelstam[s, t, u, p1, p2, -k1, -k2, 0, 0, 0, 0]; 
```

## Evaluate the amplitudes

```mathematica
$KeepLogDivergentScalelessIntegrals = True; 
```

```mathematica
ampLoop[1] = (FCTraceFactor /@ DotSimplify[#1, 
              Expanding -> False] & ) /@ Join[ampLoop[0][[1 ;; 4]], 
         Nf*ampLoop[0][[5 ;; 5]]]; 
```

```mathematica
ampTree[1] = (FCTraceFactor /@ DotSimplify[#1, 
              Expanding -> False] & ) /@ ampTree[0]; 
```

```mathematica
ampLoopCT[1] = (FCTraceFactor /@ DotSimplify[#1, 
              Expanding -> False] & ) /@ ampLoopCT[0]; 
```

```mathematica
evalFuSimple[ex_] := (Collect2[#1, {A0, B0, C0, D0}, 
            Factoring -> Function[x, Factor2[TrickMandelstam[x, 
                    {s, t, u, 0}]]]] & )[FeynAmpDenominatorExplicit[
         (#1 /. (h : A0 | B0 | C0 | D0)[x__] :> 
         TrickMandelstam[h[x], 
                    {s, t, u, 0}] & )[Contract[DiracSimplify[
               (TID[#1, q, ToPaVe -> True] & )[DiracSimplify[
                   Contract[ex]]]]]]]]; 
```

```mathematica
AbsoluteTiming[ampLoop[2] = evalFuSimple /@ ampLoop[1]; ]
```

![0hf7m556tot8b](img/0hf7m556tot8b.svg)

```mathematica
ampTree[2] = (FCCanonicalizeDummyIndices[#1, 
          LorentzIndexNames -> {mu}] & )[FeynAmpDenominatorExplicit[
       DiracSimplify[Contract[Total[ampTree[1]]]]]]
```

![11j2xyvnjbfej](img/11j2xyvnjbfej.svg)

Obtain the Born-virtual interference term

```mathematica
AbsoluteTiming[bornVirtualUnrenormalized[0] = 
       (Collect2[#1, B0, C0, D0] & )[
         (TrickMandelstam[#1, {s, t, u, 0}] & )[
           FRH[DiracSimplify[(FermionSpinSum[#1, ExtraFactor -> 
                        1/2^2] & )[Collect2[Total[ampLoop[2]], Spinor, 
                     LorentzIndex, IsolateNames -> KK]*
         ComplexConjugate[
                     ampTree[2]]]]]]]; ]
```

![0lg7j9fe2tvm7](img/0lg7j9fe2tvm7.svg)

The explicit expressions for the PaVe functions can be obtained e.g. using Package-X / PaXEvaluate

```mathematica
PaVeEvalRules = {B0[0, 0, 0] -> -(16*EpsilonIR*Pi^4)^(-1) + 
           1/(16*EpsilonUV*Pi^4), B0[s_, 0, 0] :> 
         1/(16*EpsilonUV*Pi^4) - (-2 + EulerGamma - Log[4*Pi] - 
                Log[-(ScaleMu^2/s)])/(16*Pi^4), 
       C0[0, s_, 0, 0, 0, 0] :> C0[0, 0, s, 0, 0, 0], 
       C0[0, 0, s_, 0, 0, 0] :> 1/(16*EpsilonIR^2*Pi^4*s) - 
           (EulerGamma - Log[4*Pi] - Log[-(ScaleMu^2/s)])/
             (16*EpsilonIR*Pi^4*s) - (-6*EulerGamma^2 + Pi^2 + 
                12*EulerGamma*Log[4*Pi] - 6*Log[4*Pi]^2 + 
                12*EulerGamma*Log[-(ScaleMu^2/s)] - 12*Log[4*Pi]*
                  Log[-(ScaleMu^2/s)] - 6*Log[-(ScaleMu^2/s)]^2)/
             (192*Pi^4*s), D0[0, 0, 0, 0, s_, t_, 0, 0, 0, 0] :> 
         1/(4*EpsilonIR^2*Pi^4*s*t) - (2*EulerGamma - 2*Log[4*Pi] - 
                Log[-(ScaleMu^2/s)] - Log[-(ScaleMu^2/t)])/
             (8*EpsilonIR*Pi^4*s*t) - (-3*EulerGamma^2 + 2*Pi^2 + 
                6*EulerGamma*Log[4*Pi] - 3*Log[4*Pi]^2 + 
                3*EulerGamma*Log[-(ScaleMu^2/s)] - 3*Log[4*Pi]*
                  Log[-(ScaleMu^2/s)] + 3*EulerGamma*
                  Log[-(ScaleMu^2/t)] - 3*Log[4*Pi]*
                  Log[-(ScaleMu^2/t)] - 3*Log[-(ScaleMu^2/s)]*
                  Log[-(ScaleMu^2/t)])/(24*Pi^4*s*t)}; 
```

```mathematica
bornVirtualUnrenormalized[1] = bornVirtualUnrenormalized[0] //. 
       PaVeEvalRules; 
```

Put together the counter-term contribution and the residue pole contribution

```mathematica
MSbarRC = {SMP["dZ_psi"] -> (-SMP["e"]^2/(16*Pi^2))*
           (1/EpsilonUV), SMP["dZ_A"] -> 
         (-Nf)*(SMP["e"]^2/(12*Pi^2))*(1/EpsilonUV)}; 
```

```mathematica
RuleRS = {dZe1 -> (-2^(-1))*SMP["dZ_A"], dZAA1 -> SMP["dZ_A"], 
       (dZf1 | dZf2)[__] -> SMP["dZ_psi"]}; 
```

```mathematica
legResidueContrib = 1 + (SMP["e"]^2/(4*Pi))*(1/(4*Pi))*
         (1/EpsilonIR); 
```

```mathematica
aux0 = (FCCanonicalizeDummyIndices[#1, LorentzIndexNames -> 
              {mu}] & )[DiracSimplify[Contract[
           FeynAmpDenominatorExplicit[Total[ampLoopCT[1]] /. 
                 RuleRS /. MSbarRC]]]]; 
```

```mathematica
ctContrib = Simplify[aux0/ampTree[2]]; 
```

```mathematica
fullCTAndResidue[0] = 
   (ctContrib + (4*(1/2))*(legResidueContrib - 1))*ampTree[2]
```

![1q1qgi825lz0a](img/1q1qgi825lz0a.svg)

Now get the interference of the counter term and residue contribution with the Born amplitude

```mathematica
bornCTAndResidue[0] = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     Simplify[DiracSimplify[
         (FermionSpinSum[#1, ExtraFactor -> 1/2^2] & )[
           fullCTAndResidue[0]*ComplexConjugate[ampTree[2]]]]]]
```

![01dg0nxcx5pff](img/01dg0nxcx5pff.svg)

For convenience, let us pull out an overall prefactor to get rid of ScaleMu, EulerGamma and some Pi's

```mathematica
aux1 = (#1 /. {EpsilonIR -> 1/SMP["Delta_IR"], 
              EpsilonUV -> 1/SMP["Delta_UV"]} & )[
       FCSplit[bornCTAndResidue[0], {EpsilonUV}]]; 
bornCTAndResidue[1] = (Collect2[#1, EpsilonUV, EpsilonIR] & )[
     Normal[(Series[#1, {EpsilonIR, 0, 0}] & )[
         Normal[(Series[#1, {EpsilonUV, 0, 0}] & )[
             FCShowEpsilon[FCReplaceD[
                   (1/Exp[EpsilonIR*(Log[4*Pi] - EulerGamma)])*
                     aux1[[1]], D -> 4 - 2*EpsilonIR] + 
                 
        FCReplaceD[(1/Exp[EpsilonUV*(Log[4*Pi] - EulerGamma)])*
                     aux1[[2]], D -> 4 - 2*EpsilonUV]]]]]]]
```

![070uclui9qec7](img/070uclui9qec7.svg)

```mathematica
aux2 = FCSplit[bornVirtualUnrenormalized[1], {EpsilonUV}]; 
bornVirtualUnrenormalized[2] = 
     (Collect2[#1, EpsilonUV, EpsilonIR] & )[
       (TrickMandelstam[#1, {s, t, u, 0}] & )[
         (#1 /. Log[-ScaleMu^2/(h : s | t | u)] :> 2*Log[ScaleMu] - 
                    Log[-h] & )[
     Normal[(Series[#1, {EpsilonIR, 0, 0}] & )[
               Normal[(Series[#1, {EpsilonUV, 0, 0}] & )[
                   Normal[(Collect2[#1, EpsilonUV, EpsilonIR] & )[
                       FCReplaceD[(1/ScaleMu^(2*EpsilonIR))*(1/
                                
                Exp[EpsilonIR*(Log[4*Pi] - EulerGamma)])*aux2[[
                                1]], D -> 4 - 2*EpsilonIR] + 
            FCReplaceD[
                           (1/ScaleMu^(2*EpsilonUV))*(1/
                Exp[EpsilonUV*
                                    (Log[4*Pi] - EulerGamma)])*
              aux2[[2]], 
                           D -> 4 - 2*EpsilonUV]]]]]]]]]]; 
```

Finally, we obtain the UV-finite but IR-divergent Born-virtual interference term

```mathematica
bornVirtualRenormalized[0] = 
   (Collect2[#1, EpsilonUV, EpsilonIR] & )[
     (TrickMandelstam[#1, {s, t, u, 0}] & )[
       bornVirtualUnrenormalized[2] + bornCTAndResidue[1]]]
```

![1voxdjplr08ag](img/1voxdjplr08ag.svg)

We can compare our O(eps^0) result to Eq. 2.22 in arXiv:hep-ph/0010075

```mathematica
ClearAll[LitA, LitATilde, auxBox6, Box6Eval, TriEval]; 
Li4 = PolyLog[4, #1] & ; 
ruleLit = {LitV -> Log[-s/u], LitW -> Log[-t/u], v -> s/u, 
       w -> t/u}; 
```

```mathematica
LitA = 4*GaugeXi*(1 - 2*Epsilon)*(u/s^2)*((2 - 3*Epsilon)*u^2 - 
            6*Epsilon*t*u + 3*(2 - Epsilon)*t^2)*Box6[s, t] - 
       4*(GaugeXi/(1 - 2*Epsilon))*(t/s^2)*
         ((4 - 12*Epsilon + 7*Epsilon^2)*t^2 - 
            6*Epsilon*(1 - 2*Epsilon)*t*u + 
            (4 - 10*Epsilon + 5*Epsilon^2)*u^2)*Tri[t] - 
       (8/((1 - 2*Epsilon)*(3 - 2*Epsilon)))*(1/s)*
         (2*Epsilon*(1 - Epsilon)*t*((1 - Epsilon)*t - Epsilon*u)*
              Nf - Epsilon*(3 - 2*Epsilon)*(2 - Epsilon + 2*Epsilon^2)*
              t*u + (1 - Epsilon)*(3 - 2*Epsilon)*
              (2 - (1 - GaugeXi)*Epsilon + 2*Epsilon^2)*t^2)*Tri[s]; 
```

```mathematica
auxBox6 = (1/2)*((LitV - LitW)^2 + Pi^2) + 
       2*Epsilon*(Li3[-v] - LitV*Li2[-v] - (1/3)*LitV^3 - 
            (Pi^2/2)*LitV) - 2*Epsilon^2*(Li4[-v] + LitW*Li3[-v] - 
            (1/2)*LitV^2*Li2[-v] - (1/8)*LitV^4 - (1/6)*LitV^3*
       LitW + 
            (1/4)*LitV^2*LitW^2 - (Pi^2/4)*LitV^2 - 
            (Pi^2/3)*LitV*LitW - 2*Zeta4); 
Box6Eval[s, t] = (u^(-1 - Epsilon)/(2*(1 - 2*Epsilon)))*
       (1 - (Pi^2/12)*Epsilon^2)*(auxBox6 + 
          (auxBox6 /. {LitW -> LitV, LitV -> LitW, v -> w, 
               w -> v})); 
Box6Eval[s, u] = Box6Eval[s, t] /. ruleLit /. {t -> u, u -> t}; 
TriEval[s_] := (-(-s)^(-1 - Epsilon)/Epsilon^2)*
     (1 - (Pi^2/12)*Epsilon^2 - (7/3)*Zeta[3]*Epsilon^3 - 
        (47/16)*Zeta4*Epsilon^4)
```

```mathematica
knownResult = (2/3)*(Nf/Epsilon)*8*((t^2 + u^2)/s^2 - 
              Epsilon) + ((LitA /. {Tri -> TriEval, 
                   Box6 -> Box6Eval} /. ruleLit) + 
            (LitA /. {Tri -> TriEval, Box6 -> Box6Eval} /. 
                   {GaugeXi -> -GaugeXi} /. ruleLit /. 
               {t -> u, u -> t})) /. GaugeXi -> 1; 
```

knownResult is the 1-loop result. Notice that is also an implicit overall prefactor prefLit from Eq. 2.8

```mathematica
prefLit = 32*(Pi^2/SMP["e"]^6); 
```

```mathematica
diff = (TrickMandelstam[#1, {s, t, u, 0}] & )[
     SimplifyPolyLog[PowerExpand[
         (TrickMandelstam[#1, {s, t, u, 0}] & )[
           Normal[Series[knownResult - prefLit*
                   (bornVirtualRenormalized[0] /. 
           EpsilonIR -> Epsilon), 
               {Epsilon, 0, 0}]]]]]]
```

![16f89yn3le0yv](img/16f89yn3le0yv.svg)

## Check the final results

```mathematica
FCCompareResults[0, diff, 
     Text -> {"\tCompare to arXiv:hep-ph/0010075:", "CORRECT.", 
         "WRONG!"}, Interrupt -> {Hold[Quit[1]], Automatic}]; 
Print["\tCPU Time used: ", Round[N[TimeUsed[], 4], 0.001], 
     " s."]; 
```

![0jiv5p570x7l4](img/0jiv5p570x7l4.svg)

![03md6psah2qw9](img/03md6psah2qw9.svg)