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
## Quantum Electrodynamics (QED)

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
| Renormalization ($\overline{\textrm{MS}}$) | X | [massive](FeynCalcExamples/QED/OneLoop/Renormalization), [massless](FeynCalcExamples/QED/OneLoop/RenormalizationMassless), [old](FeynCalcExamples/QED/OneLoop/Renormalization2) | [$\gamma \to \gamma$](FeynCalcExamples/QED/TwoLoops/Renormalization-GaGa), [$\ell \to \ell$](FeynCalcExamples/QED/TwoLoops/Renormalization-LeAle), [massless $\ell \to \ell$](FeynCalcExamples/QED/TwoLoops/Renormalization-LeAle-Massless) |
| $e^+ e^- \to e^+ e^-$ | [link](FeynCalcExamples/QED/Tree/ElAel-ElAel) | [link](FeynCalcExamples/QED/OneLoop/ElAel-ElAel) | |
| $e^+ e^- \to \mu^+ \mu^-$ | [link](FeynCalcExamples/QED/Tree/ElAel-MuAmu) | [link](FeynCalcExamples/QED/OneLoop/ElAel-MuAmu) | |
| $e^+ e^- \to \gamma \gamma$ | [link](FeynCalcExamples/QED/Tree/ElAel-GaGa) | | |
| $e^- \mu^+ \to e^- \mu^+$ | [link](FeynCalcExamples/QED/Tree/ElMu-ElMu) | | |
| $e^- \gamma \to e^- \gamma$ | [link](FeynCalcExamples/QED/Tree/ElGa-ElGa) | | |
| $e^- e^- \to e^- e^-$ | [link](FeynCalcExamples/QED/Tree/ElEl-ElEl) | | |
| $e^- \to \gamma e^-$ | X | [$g-2$](FeynCalcExamples/QED/OneLoop/El-GaEl) | |
|$\gamma$ | X | [Furry's theorem](FeynCalcExamples/QED/OneLoop/Ga) | |
|$\gamma \to \gamma \gamma$ | X | [Furry's theorem](FeynCalcExamples/QED/OneLoop/Ga-GaGa) | |
|$\gamma \to \gamma \gamma \gamma \gamma$ | X | [Furry's theorem](FeynCalcExamples/QED/OneLoop/Ga-GaGaGaGa) | |
|$\gamma \to \mu^+ \mu^-$ | [link](FeynCalcExamples/QED/Tree/Ga-MuAmu) |  | |
|$\pi \to \gamma \gamma$ | X | [Axial anomaly](FeynCalcExamples/QED/OneLoop/Pi-GaGa) | |
|:---:|:---:|:---:|:---:|


## Quantum chromodynamics (QCD)


|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
| Renormalization ($\overline{\textrm{MS}}$) | X | [massive](FeynCalcExamples/QCD/OneLoop/Renormalization), [massless](FeynCalcExamples/QCD/OneLoop/RenormalizationMassless), [old massive](FeynCalcExamples/QCD/OneLoop/Renormalization2), [old massless](FeynCalcExamples/QCD/OneLoop/RenormalizationMassless2) | [$g \to g$](FeynCalcExamples/QCD/TwoLoops/Renormalization-GlGl), [$u_g \to u_g$](FeynCalcExamples/QCD/TwoLoops/Renormalization-GhGh), [$q \to q$](FeynCalcExamples/QCD/TwoLoops/Renormalization-QQbar), [massless $q \to q$](FeynCalcExamples/QCD/TwoLoops/Renormalization-QQbar-Massless), [$q g \to q$](FeynCalcExamples/QCD/TwoLoops/Renormalization-QQbarGl) |
| $g g \to H$| X | | [link](FeynCalcExamples/QCD/TwoLoops/GlGl-H) |
| $g g \to H H$| X | [link](FeynCalcExamples/QCD/OneLoop/GlGl-HH) |  |
| $H \to g g$| X | [link](FeynCalcExamples/QCD/OneLoop/H-GlGl) |  |
| $g \to g$| X | [BG field gauge](FeynCalcExamples/QCD/OneLoop/Gl-Gl-BackgroundFieldGauge) | [link](FeynCalcExamples/QCD/TwoLoops/Gl-Gl) |
| $u_g \to u_g$| X |  | [link](FeynCalcExamples/QCD/TwoLoops/Gh-Gh) |
| $g \to g g$| [link](FeynCalcExamples/QCD/OneLoop/Gl-GlGl) | | |
| massless NS Adler function | X | [link](FeynCalcExamples/QCD/OneLoop/Adler-massless) | [link](FeynCalcExamples/QCD/TwoLoops/Adler-massless) |
| $e^+ e^- \to  q \bar{q}$| [link](FeynCalcExamples/QCD/Tree/ElAel-QQbar) | | |
| $\gamma g \to  q \bar{q}$| [link](FeynCalcExamples/QCD/Tree/GaGl-QQbar) | | |
| $\gamma \to  q \bar{q}$| [link](FeynCalcExamples/QCD/Tree/Ga-QQbar) | | |
| $\gamma \to  q \bar{q} g$| [link](FeynCalcExamples/QCD/Tree/Ga-QQbarGl) | | |
| $g g \to g g$| [link](FeynCalcExamples/QCD/Tree/GlGl-GlGl) | | |
| $g g\to  q \bar{q}$| [link](FeynCalcExamples/QCD/Tree/GlGl-QQbar) | | |
| $\mu^+ \mu^- \to  q \bar{q}$| [link](FeynCalcExamples/QCD/Tree/MuAmu-QQbar) | | |
| $q \gamma \to g q$| [link](FeynCalcExamples/QCD/Tree/QGa-GlQ) | | |
| $q g \to q g$| [v1](FeynCalcExamples/QCD/Tree/QGl-QGl), [v2](FeynCalcExamples/QCD/Tree/QGl-QGl-2) | | |
| $q_i \bar{q}_i \to q_i \bar{q}_i$| [link](FeynCalcExamples/QCD/Tree/QiQibar-QiQibar) | | |
| $q_i \bar{q}_i \to q_j \bar{q}_j$| [link](FeynCalcExamples/QCD/Tree/QiQibar-QjQjbar) |  |
| $q_i q_i \to q_i q_i$| [link](FeynCalcExamples/QCD/Tree/QiQi-QiQi) |  |
| $q_i \bar{q}_j \to q_i \bar{q}_j$| [link](FeynCalcExamples/QCD/Tree/QiQjbar-QiQjbar) |  |
| $q_i q_j \to q_i q_j$| [link](FeynCalcExamples/QCD/Tree/QiQj-QiQj) |  |
| $q \bar{q} \to e^+ e^-$| [link](FeynCalcExamples/QCD/Tree/QQbar-ElAel) |  |
| $q \bar{q} \to \gamma \gamma$| [link](FeynCalcExamples/QCD/Tree/QQbar-GaGa) |  |
| $q \bar{q} \to \gamma g$| [link](FeynCalcExamples/QCD/Tree/QQbar-GaGl) |  |
| $q \bar{q} \to g g$| [v1](FeynCalcExamples/QCD/Tree/QQbar-GlGl), [v2](FeynCalcExamples/QCD/Tree/QQbar-GlGl-2) |  |
| $q \bar{q} \to \mu^+ \mu^-$| [link](FeynCalcExamples/QCD/Tree/QQbar-MuAmu) |   |
| Soft function $\gamma \to q\bar{q}$ | [link](FeynCalcExamples/QCD/Tree/Ga-QQbar-SoftFunction) |   |
|:---:|:---:|:---:|:---:|

## Electroweak theory (EW)

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
| $q \bar{q} \to Z^0 Z^0$ | [link](FeynCalcExamples/EW/Tree/QQbar-ZZ) |   | |
| $H \to f \bar{f}$ | [link](FeynCalcExamples/EW/Tree/H-FFbar) |   | |
| $H \to W^+ W^-$ | [link](FeynCalcExamples/EW/Tree/H-WW) |   | |
| $H \to Z^0 Z^0$ | [link](FeynCalcExamples/EW/Tree/H-ZZ) |   | |
| $W \to q_i \bar{q}_j$ | [link](FeynCalcExamples/EW/Tree/W-QiQjbar) | | |
| $Z \to f \bar{f}$ | [link](FeynCalcExamples/EW/Tree/Z-FFbar) | | |
| ${\bar{\nu}}_e e^- \to {\bar{\nu}}_{\mu} \mu$ | [link](FeynCalcExamples/EW/Tree/AnelEl-AnmuMu) | |  |
| $\bar{\nu}_e e^- \to \bar{q}_u q_d$ | [link](FeynCalcExamples/EW/Tree/AnelEl-QubarQd) | | |
| $\bar{\nu}_e e^- \to W^+ W^-$ | [link](FeynCalcExamples/EW/Tree/AnelEl-WW) |   | |
| $\bar{\nu}_e e^- \to Z^0 Z^0$ | [link](FeynCalcExamples/EW/Tree/AnelEl-ZZ) |   | |
| $e^- \nu_{\mu} \to \mu \nu_{e}$ | [link](FeynCalcExamples/EW/Tree/ElNmu-MuNel) |   | |
| $\nu_l q_d \to l q_u$ | [link](FeynCalcExamples/EW/Tree/NleQdt-LeQut) |   | |
| $\bar{q}_d q_u \to \bar{\nu}_e \nu_e$ | [link](FeynCalcExamples/EW/Tree/QuQdbar-AelNel) |   | |
| $\bar{q}_u q_d \to \bar{\nu}_e \nu_e$ | [link](FeynCalcExamples/EW/Tree/QutbarQdt-NelAnel) |   | |
| $W^- \to e^- \bar{\nu}_e$  | [link](FeynCalcExamples/EW/Tree/W-ElAnel) |   | |
| $$\mu \to e^- {\bar{\nu}}_{e} \nu_{\mu}$$ | [link](FeynCalcExamples/EW/Tree/Mu-ElAnelNmu) |   | |
| $t \to W^+ b$ | [link](FeynCalcExamples/EW/Tree/Qt-QbW) |   | |
|:---:|:---:|:---:|:---:|



## $\phi^3$ theory

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
 Renormalization ($\overline{\textrm{MS}}$) | X | [link](FeynCalcExamples/Phi3/OneLoop/Renormalization) |   | |
|:---:|:---:|:---:|:---:|

## $\phi^4$ theory

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
 Renormalization ($\overline{\textrm{MS}}$) | X | [link](FeynCalcExamples/Phi4/OneLoop/Renormalization) | [$\phi \to \phi$](FeynCalcExamples/Phi4/TwoLoops/Renormalization-SS), [$\phi \phi \to \phi \phi$](FeynCalcExamples/Phi4/TwoLoops/Renormalization-SSSS)  | |
| $\phi \phi \to \phi \phi$ |  | [link](FeynCalcExamples/Phi4/OneLoop/PhiPhi-PhiPhi) | |
|:---:|:---:|:---:|:---:|

## Scalar Yukawa theory

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
 Renormalization ($\overline{\textrm{MS}}$) | X | [link](FeynCalcExamples/YukawaS/OneLoop/Renormalization) | | |
|:---:|:---:|:---:|:---:|

## Pseudoscalar Yukawa theory

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
 Renormalization ($\overline{\textrm{MS}}$) | X | [link](FeynCalcExamples/YukawaPS/OneLoop/Renormalization) | | |
|:---:|:---:|:---:|:---:|

## Minimal Supersymmetric Standard Model (MSSM)

|Process |Tree level |1-loop level|2-loop level|
|:---:|:---:|:---:|:---:|
| $\tilde{\nu}_e e^- \to \tilde{\nu}_e e^-$| [link](FeynCalcExamples/MSSM/Tree/MnelEl-MnelEl) |  | |
|:---:|:---:|:---:|:---:|

## Topology identification and minimization

|Process |1-loop level|2-loop level|3-loop level|
|:---:|:---:|:---:|:---:|
| $B_c \to \eta_c$| X | [link](FeynCalcExamples/TopologyIdentification/TwoLoops/B-EtaC) |  |
|:---:|:---:|:---:|:---:|

## Acknowledgement

To generate Markdown files out of FeynCalc notebooks we employed the [M2MD](https://github.com/kubaPod/M2MD)
package by [Kuba Podkalicki](https://github.com/kubaPod). We thank the author for helpful advice on using the package in a proper way.

