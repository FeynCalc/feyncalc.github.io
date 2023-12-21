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
|:-------------:|:-------------:|
| | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamplesMD/Phi3/OneLoop/Renormalization) |
|:-------------:|:-------------:|


## $\phi^4$ theory

|Tree level |1-loop level|
|:-------------:|:-------------:|
| | [$\phi \phi \to \phi \phi$](FeynCalcExamplesMD/Phi4/OneLoop/PhiPhi-PhiPhi) |
| | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamplesMD/Phi4/OneLoop/Renormalization) |
|:-------------:|:-------------:|



## Quantum electrodynamics (QED)

|Tree level |1-loop level|2-loop level|
|:-------------:|:-------------:|
| [$e^+ e^- \to e^+ e^-$](FeynCalcExamplesMD/QED/Tree/ElAel-ElAel) | [$e^+ e^- \to e^+ e^-$](FeynCalcExamplesMD/QED/OneLoop/ElAel-ElAel) | [$e^- \to e^-$](FeynCalcExamplesMD/QED/TwoLoops/El-El) |
| [$e^+ e^- \to \gamma \gamma$](FeynCalcExamplesMD/QED/Tree/ElAel-GaGa) | [$e^- \to e^-$](FeynCalcExamplesMD/QED/OneLoop/El-El) | [$\gamma \to \gamma$](FeynCalcExamplesMD/QED/TwoLoops/Ga-Ga) |
| [$e^+ e^- \to \mu^+ \mu^-$](FeynCalcExamplesMD/QED/Tree/ElAel-MuAmu) | [$e^- \to \gamma e^-$](FeynCalcExamplesMD/QED/OneLoop/El-GaEl) | |
| [$e^- e^- \to e^- e^-$](FeynCalcExamplesMD/QED/Tree/ElEl-ElEl) | [$\gamma$](FeynCalcExamplesMD/QED/OneLoop/Ga) | |
| [$e^- \gamma \to e^- \gamma$](FeynCalcExamplesMD/QED/Tree/ElGa-ElGa) | [$\gamma \to \gamma$](FeynCalcExamplesMD/QED/OneLoop/Ga-Ga) | |
| [$e^- \mu^+ \to e^- \mu^+$](FeynCalcExamplesMD/QED/Tree/ElMu-ElMu) | [$\gamma \to \gamma \gamma$](FeynCalcExamplesMD/QED/OneLoop/Ga-GaGa) | |
| [$\gamma \to \mu^+ \mu^-$](FeynCalcExamplesMD/QED/Tree/Ga-MuAmu) | [$\gamma \to \gamma \gamma \gamma \gamma$](FeynCalcExamplesMD/QED/OneLoop/Ga-GaGaGaGa) | |
| | [$e^+ e^- \to \mu^+ \mu^-$](FeynCalcExamplesMD/QED/OneLoop/ElAel-MuAmu) | |
| | [$\pi \to \gamma \gamma$](FeynCalcExamplesMD/QED/OneLoop/Pi-GaGa) | |
| | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamplesMD/QED/OneLoop/Renormalization) | |
|:-------------:|:-------------:|

## Quantum chromodynamics (QCD)

|Tree level |1-loop level| 2-loop level|
|:-------------:|:-------------:|
| [$e^+ e^- \to  q \bar{q}$](FeynCalcExamplesMD/QCD/Tree/ElAel-QQbar) | [$u_g \to u_g$](FeynCalcExamplesMD/QCD/OneLoop/Gh-Gh) | [$u_g g \to u_g$](FeynCalcExamplesMD/QCD/TwoLoops/Gh-Gh) |  
| [$\gamma g \to  q \bar{q}$](FeynCalcExamplesMD/QCD/Tree/GaGl-QQbar) | [$u_g g \to u_g$](FeynCalcExamplesMD/QCD/OneLoop/GhGl-Gh) | [$g \to g$](FeynCalcExamplesMD/QCD/TwoLoops/Gl-Gl) | 
| [$\gamma \to  q \bar{q}$](FeynCalcExamplesMD/QCD/Tree/Ga-QQbar) | [$g \to g$](FeynCalcExamplesMD/QCD/OneLoop/Gl-Gl) |  | 
| [$\gamma \to  q \bar{q} g$](FeynCalcExamplesMD/QCD/Tree/Ga-QQbarGl) | [$g \to g$ (background field gauge)](FeynCalcExamplesMD/QCD/OneLoop/Gl-Gl-BackgroundFieldGauge) |  |
| [$g g \to g g$](FeynCalcExamplesMD/QCD/Tree/GlGl-GlGl) | [$g \to g g$](FeynCalcExamplesMD/QCD/OneLoop/Gl-GlGl) |  |
| [$g g\to  q \bar{q}$](FeynCalcExamplesMD/QCD/Tree/GlGl-QQbar) | [$q \to q$](FeynCalcExamplesMD/QCD/OneLoop/Q-Q) |  |
| [$\mu^+ \mu^- \to  q \bar{q}$](FeynCalcExamplesMD/QCD/Tree/MuAmu-QQbar) | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$)](FeynCalcExamplesMD/QCD/OneLoop/Renormalization) |  |
| [$q \gamma \to g q$](FeynCalcExamplesMD/QCD/Tree/QGa-GlQ) | [Renormalization ($\textrm{MS}, \overline{\textrm{MS}}$), massless quarks](FeynCalcExamplesMD/QCD/OneLoop/RenormalizationMassless) |  |
| [$q g \to q g$](FeynCalcExamplesMD/QCD/Tree/QGl-QGl) | |
| [$q_i \bar{q}_i \to q_i \bar{q}_i$](FeynCalcExamplesMD/QCD/Tree/QiQibar-QiQibar) |  |
| [$q_i \bar{q}_i \to q_j \bar{q}_j$](FeynCalcExamplesMD/QCD/Tree/QiQibar-QjQjbar) |  |
| [$q_i q_i \to q_i q_i$](FeynCalcExamplesMD/QCD/Tree/QiQi-QiQi) |  |
| [$q_i \bar{q}_j \to q_i \bar{q}_j$](FeynCalcExamplesMD/QCD/Tree/QiQjbar-QiQjbar) |  |
| [$q_i q_j \to q_i q_j$](FeynCalcExamplesMD/QCD/Tree/QiQj-QiQj) |  |
| [$q \bar{q} \to e^+ e^-$](FeynCalcExamplesMD/QCD/Tree/QQbar-ElAel) |  |
| [$q \bar{q} \to \gamma \gamma$](FeynCalcExamplesMD/QCD/Tree/QQbar-GaGa) |  |
| [$q \bar{q} \to \gamma g$](FeynCalcExamplesMD/QCD/Tree/QQbar-GaGl) |  |
| [$q \bar{q} \to g g$](FeynCalcExamplesMD/QCD/Tree/QQbar-GlGl) |  |
| [$q \bar{q} \to \mu^+ \mu^-$](FeynCalcExamplesMD/QCD/Tree/QQbar-MuAmu) |   |
|:-------------:|:-------------:|

## Electroweak theory (EW)

|Tree level |1-loop level|
|:-------------:|:-------------:|
| [$$\bar{\nu}_e e^- \to \bar{\nu}_{\mu} \mu$$](FeynCalcExamplesMD/EW/Tree/AnelEl-AnmuMu) |  [$H \to g g$](FeynCalcExamplesMD/EW/OneLoop/H-GG) | 
| [$\bar{\nu}_e e^- \to \bar{q}_u q_d$](FeynCalcExamplesMD/EW/Tree/AnelEl-QubarQd) |   |
| [$\bar{\nu}_e e^- \to W^+ W^-$](FeynCalcExamplesMD/EW/Tree/AnelEl-WW) |   |
| [$\bar{\nu}_e e^- \to Z^0 Z^0$](FeynCalcExamplesMD/EW/Tree/AnelEl-ZZ) |   |
| [$e^- \nu_{\mu} \to \mu \nu_{e}$](FeynCalcExamplesMD/EW/Tree/ElNmu-MuNel) |   |
| [$H \to f \bar{f}$](FeynCalcExamplesMD/EW/Tree/H-FFbar) |   |
| [$H \to W^+ W^-$](FeynCalcExamplesMD/EW/Tree/H-WW) |   |
| [$H \to Z^0 Z^0$](FeynCalcExamplesMD/EW/Tree/H-ZZ) |   |
| [$$\mu \to e^- \bar{\nu}_{e} \nu_{\mu}$$](FeynCalcExamplesMD/EW/Tree/Mu-ElAnelNmu) |   |
| [$\nu_l q_d \to l q_u$](FeynCalcExamplesMD/EW/Tree/NleQdt-LeQut) |   |
| [$q \bar{q} \to Z^0 Z^0$](FeynCalcExamplesMD/EW/Tree/QQbar-ZZ) |   |
| [$t \to W^+ b$](FeynCalcExamplesMD/EW/Tree/Qt-QbW) |   |
| [$\bar{q}_d q_u \to \bar{\nu}_e \nu_e$](FeynCalcExamplesMD/EW/Tree/QuQdbar-AelNel) |   |
| [$\bar{q}_u q_d \to \bar{\nu}_e \nu_e$](FeynCalcExamplesMD/EW/Tree/QutbarQdt-NelAnel) |   |
| [$W^- \to e^- \bar{\nu}_e $](FeynCalcExamplesMD/EW/Tree/W-ElAnel) |   |
| [$W \to q_i \bar{q}_j$](FeynCalcExamplesMD/EW/Tree/W-QiQjbar) |   |
| [$Z \to f \bar{f}$](FeynCalcExamplesMD/EW/Tree/Z-FFbar) |   |
|:-------------:|:-------------:|

## Minimal Supersymmetric Standard Model (MSSM)

|Tree level |1-loop level|
|:-------------:|:-------------:|
| [$$\tilde{\nu}_e e^- \to \tilde{\nu}_e e^-$$](FeynCalcExamplesMD/MSSM/Tree/MnelEl-MnelEl) |  |
|:-------------:|:-------------:|

<a name="ack"/>
## Acknowledgement

To generate Markdown files out of FeynCalc notebooks we employed the [M2MD](https://github.com/kubaPod/M2MD)
package by [Kuba Podkalicki](https://github.com/kubaPod). We thank the author for helpful advice on using the package in a proper way.
Our scripts for running M2MD can be found [here](https://github.com/FeynCalc/feyncalc/tree/master/FeynCalc/Examples/Export).






