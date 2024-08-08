---
layout: default
title: Example gallery
permalink: /examples
---

Here we collect examples of interesting calculations from textbooks and papers that can be reproduced using
_FeynCalc_. The HTML pages below were [generated](#ack) from [example notebooks](https://github.com/FeynCalc/feyncalc/tree/hotfix-stable/FeynCalc/Examples) shipped with the current stable version
of Feyncalc. If you have FeynCalc installed, these files are located in 

```mathematica 
FileNameJoin[{$UserBaseDirectory, "Applications", "FeynCalc", "Examples"}]
```

## $\phi^3$ theory

|Tree level |1-loop level|
|:---:|:---:|
|  | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamples/Phi3/OneLoop/Renormalization) |
|:---:|:---:|


## $\phi^4$ theory

|Tree level |1-loop level|
|:---:|:---:|
|  | [$\phi \phi \to \phi \phi$](FeynCalcExamples/Phi4/OneLoop/PhiPhi-PhiPhi) |
|  | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamples/Phi4/OneLoop/Renormalization) |
|:---:|:---:|



## Quantum Electrodynamics (QED)

|Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|
| [$e^+ e^- \to e^+ e^-$](FeynCalcExamples/QED/Tree/ElAel-ElAel) | [$e^+ e^- \to e^+ e^-$ v.1](FeynCalcExamples/QED/OneLoop/ElAel-ElAel) | [$e^- \to e^-$](FeynCalcExamples/QED/TwoLoops/El-El) |
| [$e^+ e^- \to \gamma \gamma$](FeynCalcExamples/QED/Tree/ElAel-GaGa) | [$e^- \to e^-$](FeynCalcExamples/QED/OneLoop/El-El) | [$\gamma \to \gamma$](FeynCalcExamples/QED/TwoLoops/Ga-Ga) |
| [$e^+ e^- \to \mu^+ \mu^-$](FeynCalcExamples/QED/Tree/ElAel-MuAmu) | [$e^- \to \gamma e^-$](FeynCalcExamples/QED/OneLoop/El-GaEl) | |
| [$e^- e^- \to e^- e^-$](FeynCalcExamples/QED/Tree/ElEl-ElEl) | [$\gamma$](FeynCalcExamples/QED/OneLoop/Ga) | |
| [$e^- \gamma \to e^- \gamma$](FeynCalcExamples/QED/Tree/ElGa-ElGa) | [$\gamma \to \gamma$](FeynCalcExamples/QED/OneLoop/Ga-Ga) | |
| [$e^- \mu^+ \to e^- \mu^+$](FeynCalcExamples/QED/Tree/ElMu-ElMu) | [$\gamma \to \gamma \gamma$](FeynCalcExamples/QED/OneLoop/Ga-GaGa) | |
| [$\gamma \to \mu^+ \mu^-$](FeynCalcExamples/QED/Tree/Ga-MuAmu) | [$\gamma \to \gamma \gamma \gamma \gamma$](FeynCalcExamples/QED/OneLoop/Ga-GaGaGaGa) | |
| | [$e^+ e^- \to \mu^+ \mu^-$ v.1](FeynCalcExamples/QED/OneLoop/ElAel-MuAmu) |  |
| | [$\pi \to \gamma \gamma$](FeynCalcExamples/QED/OneLoop/Pi-GaGa) | |
| | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamples/QED/OneLoop/Renormalization) | |
| | [$e^+ e^- \to e^+ e^-$ v.2](FeynCalcExamples/QED/OneLoop/ElAel-ElAel2) | |
| | [$e^+ e^- \to \mu^+ \mu^-$ v.2](FeynCalcExamples/QED/OneLoop/ElAel-MuAmu2) | |
| | [$\pi \to \gamma \gamma$](FeynCalcExamples/QED/OneLoop/PiToGaGa) | |
|:---:|:---:|:---:|

## Quantum chromodynamics (QCD)

|Tree level |1-loop level| 2-loop level|
|:---:|:---:|:---:|
| [$e^+ e^- \to  q \bar{q}$](FeynCalcExamples/QCD/Tree/ElAel-QQbar) | [$u_g \to u_g$](FeynCalcExamples/QCD/OneLoop/Gh-Gh) | [$u_g \to u_g$ v.1 ](FeynCalcExamples/QCD/TwoLoops/Gh-Gh) |  
| [$\gamma g \to  q \bar{q}$](FeynCalcExamples/QCD/Tree/GaGl-QQbar) | [$u_g g \to u_g$](FeynCalcExamples/QCD/OneLoop/GhGl-Gh) |  [$u_g \to u_g$ v.2](FeynCalcExamples/QCD/TwoLoops/Gh-Gh-2) | 
| [$\gamma \to  q \bar{q}$](FeynCalcExamples/QCD/Tree/Ga-QQbar) | [$g \to g$](FeynCalcExamples/QCD/OneLoop/Gl-Gl) | [$g \to g$](FeynCalcExamples/QCD/TwoLoops/Gl-Gl) | 
| [$\gamma \to  q \bar{q} g$](FeynCalcExamples/QCD/Tree/Ga-QQbarGl) | [$g \to g$ (background field gauge)](FeynCalcExamples/QCD/OneLoop/Gl-Gl-BackgroundFieldGauge) |  |
| [$g g \to g g$](FeynCalcExamples/QCD/Tree/GlGl-GlGl) | [$g \to g g$](FeynCalcExamples/QCD/OneLoop/Gl-GlGl) |  |
| [$g g\to  q \bar{q}$](FeynCalcExamples/QCD/Tree/GlGl-QQbar) | [$q \to q$](FeynCalcExamples/QCD/OneLoop/Q-Q) |  |
| [$\mu^+ \mu^- \to  q \bar{q}$](FeynCalcExamples/QCD/Tree/MuAmu-QQbar) | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamples/QCD/OneLoop/Renormalization) |  |
| [$q \gamma \to g q$](FeynCalcExamples/QCD/Tree/QGa-GlQ) | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$), massless quarks](FeynCalcExamples/QCD/OneLoop/RenormalizationMassless) |  |
| [$q g \to q g$ v.1](FeynCalcExamples/QCD/Tree/QGl-QGl) | |
| [$q_i \bar{q}_i \to q_i \bar{q}_i$](FeynCalcExamples/QCD/Tree/QiQibar-QiQibar) |  |
| [$q_i \bar{q}_i \to q_j \bar{q}_j$](FeynCalcExamples/QCD/Tree/QiQibar-QjQjbar) |  |
| [$q_i q_i \to q_i q_i$](FeynCalcExamples/QCD/Tree/QiQi-QiQi) |  |
| [$q_i \bar{q}_j \to q_i \bar{q}_j$](FeynCalcExamples/QCD/Tree/QiQjbar-QiQjbar) |  |
| [$q_i q_j \to q_i q_j$](FeynCalcExamples/QCD/Tree/QiQj-QiQj) |  |
| [$q \bar{q} \to e^+ e^-$](FeynCalcExamples/QCD/Tree/QQbar-ElAel) |  |
| [$q \bar{q} \to \gamma \gamma$](FeynCalcExamples/QCD/Tree/QQbar-GaGa) |  |
| [$q \bar{q} \to \gamma g$](FeynCalcExamples/QCD/Tree/QQbar-GaGl) |  |
| [$q \bar{q} \to g g$](FeynCalcExamples/QCD/Tree/QQbar-GlGl) |  |
| [$q \bar{q} \to \mu^+ \mu^-$](FeynCalcExamples/QCD/Tree/QQbar-MuAmu) |   |
| [ Soft function $\gamma \to q \bar{q} $](FeynCalcExamples/QCD/Tree/Ga-QQbar-SoftFunction) |   |
| [$q g \to q g$ v.2](FeynCalcExamples/QCD/Tree/QGl-QGl-2) | |
| [$q \bar{q} \to g g$ v.2 ](FeynCalcExamples/QCD/Tree/QQbar-GlGl-2) |  |
|:---:|:---:|:---:|

## Electroweak theory (EW)

|Tree level |1-loop level|
|:---:|:---:
| [$$\bar{\nu}_e e^- \to \bar{\nu}_{\mu} \mu$$](FeynCalcExamples/EW/Tree/AnelEl-AnmuMu) |  [$H \to g g$](FeynCalcExamples/EW/OneLoop/H-GG) | 
| [$\bar{\nu}_e e^- \to \bar{q}_u q_d$](FeynCalcExamples/EW/Tree/AnelEl-QubarQd) |  [$\Pi \to \gamma \gamma$](FeynCalcExamples/QED/OneLoop/Pi-GaGa)  |
| [$\bar{\nu}_e e^- \to W^+ W^-$](FeynCalcExamples/EW/Tree/AnelEl-WW) |   |
| [$\bar{\nu}_e e^- \to Z^0 Z^0$](FeynCalcExamples/EW/Tree/AnelEl-ZZ) |   |
| [$e^- \nu_{\mu} \to \mu \nu_{e}$](FeynCalcExamples/EW/Tree/ElNmu-MuNel) |   |
| [$H \to f \bar{f}$](FeynCalcExamples/EW/Tree/H-FFbar) |   |
| [$H \to W^+ W^-$](FeynCalcExamples/EW/Tree/H-WW) |   |
| [$H \to Z^0 Z^0$](FeynCalcExamples/EW/Tree/H-ZZ) |   |
| [$$\mu \to e^- \bar{\nu}_{e} \nu_{\mu}$$](FeynCalcExamples/EW/Tree/Mu-ElAnelNmu) |   |
| [$\nu_l q_d \to l q_u$](FeynCalcExamples/EW/Tree/NleQdt-LeQut) |   |
| [$q \bar{q} \to Z^0 Z^0$](FeynCalcExamples/EW/Tree/QQbar-ZZ) |   |
| [$t \to W^+ b$](FeynCalcExamples/EW/Tree/Qt-QbW) |   |
| [$\bar{q}_d q_u \to \bar{\nu}_e \nu_e$](FeynCalcExamples/EW/Tree/QuQdbar-AelNel) |   |
| [$\bar{q}_u q_d \to \bar{\nu}_e \nu_e$](FeynCalcExamples/EW/Tree/QutbarQdt-NelAnel) |   |
| [$W^- \to e^- \bar{\nu}_e $](FeynCalcExamples/EW/Tree/W-ElAnel) |   |
| [$W \to q_i \bar{q}_j$](FeynCalcExamples/EW/Tree/W-QiQjbar) |   |
| [$Z \to f \bar{f}$](FeynCalcExamples/EW/Tree/Z-FFbar) |   |
|:---:|:---:|

## Yukawa theory

|Tree level |1-loop level|
|:---:|:---:|
|  | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamples/Yukawa/OneLoop/Renormalization) |
|:---:|:---:|


## Minimal Supersymmetric Standard Model (MSSM)

|Tree level |1-loop level|
|:---:|:---:|
| [$$\tilde{\nu}_e e^- \to \tilde{\nu}_e e^-$$](FeynCalcExamples/MSSM/Tree/MnelEl-MnelEl) |  |
|:---:|:---:|

## Topology identification and minimization

|1-loop level |2-loop level| 3-loop level|
|:---:|:---:|:---:|
| | [$$ B_c \to \eta_c$$](FeynCalcExamples/TopologyIdentification/TwoLoops/B-EtaC) |  |
|:---:|:---:|:---:|

## FeynRules models for FeynArts

|Models |
|:---:|
| [$$\phi^3 $$](FeynCalcExamples/FeynRules/Phi3/GenerateModelPhi3) |
| [$$\phi^4 $$](FeynCalcExamples/FeynRules/Phi4/GenerateModelPhi4) |
| [$$\phi^3 + \phi^4 $$](FeynCalcExamples/FeynRules/Phi34/GenerateModelPhi34) |
| [Euler-Heisenberg EFT](FeynCalcExamples/FeynRules/EulerHeisenberg/GenerateModelEulerHeisenberg) |
| [QED](FeynCalcExamples/FeynRules/QED/GenerateModelQED) |
| [QCD](FeynCalcExamples/FeynRules/QCD/GenerateModelQCD) |
| [QCD in background field gauge](FeynCalcExamples/FeynRules/QCDBGF/GenerateModelQCDBGF) |
| [Standard Model](FeynCalcExamples/FeynRules/SM/GenerateModelSM) |
| [Yukawa](FeynCalcExamples/FeynRules/Yukawa/GenerateModelYukawa) |
|:---:|

<a name="ack"/>
## Acknowledgement

To generate Markdown files out of FeynCalc notebooks we employed the [M2MD](https://github.com/kubaPod/M2MD)
package by [Kuba Podkalicki](https://github.com/kubaPod). We thank the author for helpful advice on using the package in a proper way.






