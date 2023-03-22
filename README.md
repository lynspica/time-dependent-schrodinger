# time-dependent-schrodinger
## **Modelling of Time-dependent Schrödinger equation using finite difference method**

This project presents the numerical analysis of time-dependent Schrödinger equation (TDSE) in 1D and 2D. Potential term is considered time-independent, to model some well-known examples of harmonic oscillator / finite-well / tunnelling / free-particle. Moreover, single-slit and double-slit experiments were numerically simulated in 2D.

The numerical analysis of 1D TDSE is currently limited with explicit methods: Explicit Euler and Crank-Nicolson methods. On the other hand, 2D TDSE is experimented with both explicit and implicit methods: Explicit Euler and Alternating Direction Implicit (ADI) methods. While the animations given below demonstrate some exemplary cases, the derivation of numerical models are reported, (see [report](https://github.com/lynspica/time-dependent-schrodinger/blob/main/report.pdf)).

The animation below shows the tunneling phenomena in 1D, where the particle is confined within 1D infinite potential well while a finite potential barrier (width = 0.5, at $x = 1$, V(x) = 0.5) is introduced.

![](https://github.com/lynspica/time-dependent-schrodinger/blob/main/figs/1D_tunneling.gif)

The animation below shows the single and double slit experiments, where the particle is confined within 2D infinite potential well while finite potential barriers creating the slit geometries (width = 0.5, slit size = 1, at $x = 1$, V(x,y) = 100) are introduced.

![](https://github.com/lynspica/time-dependent-schrodinger/blob/main/figs/2D_slits.gif)
