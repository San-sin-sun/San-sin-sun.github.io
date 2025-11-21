---
title: Biological Physics
date: 2024-10-28 16:26:57
categories:
  - 课程笔记
tags:
math: true
---

## Intro-Biomolecules

## Statistical-Physics-Basics

### Basics

**Microstate**:
M particles, N levels: $\Omega = M^N$
**Macrostate**:
Histogram of the number of particles per energy level
Macrostate{3, 1, 5}:

- 3 particles are in level 1
- 1 particle is in level 2
- 5 particles are in level 3

Fundamental concept of statistical Physics

- Describe a system only as a MACROSTATE = Histogram
- It is not relevant which individual particle is in what level
- Compared to the exact description of a system (MICROSTATE), information is lost. A histogram contains LESS information than the data set that generated it 
- Cannot construct the exact(MICROSTATE) of a system from the MACROSTATE because the MACROSTATE is just a distribution
- The amount of the missing information is measureed by the statistical weight $\Omega$ of the Macrostate, or $\ln(\Omega)$

**Stastical weight of a macrostate $\Omega$**

How many microstates belong to the macrostate
eg: {2,0,0,1,0} has 3 microstates

$\Omega = N!$(所有粒子能量均不相同)

***Degeneracy***

${1,0,0,2}: \Omega = \frac{3!}{2!}$

### Entropy

$S = k_B \ln \Omega, \Omega = \frac{N!}{n_0!n_1!\cdots}$

**Stirling's theorem**
$$N! \approx \frac{N^N}{e^N}\sqrt{2\pi N}$$
$$\ln N! \approx N \ln N -N$$

Entropy S is a property of a macrostate of a system, gives the level of disorder of that macrostate

## Boltzman-Free-Energy

### Boltzman

$\ln \Omega = \ln \frac{N!}{\prod_{i}n_i!} = N\ln N - \sum_{i}(n_i\ln n_i)$
Maximum of entropy occurs for $dS = d(k_B \ln \Omega)=0$
$$d\ln\Omega=-\sum_{i}\ln n_i dn_i=0$$
$\sum_i n_i = N, \sum_i dn_i=0$
Use Lagrange's method of undertimined multipliers: $n_i = \frac{N}{M}$

#### Microcanonical ensemble

$$\sum_i \epsilon_in_i = E$$
$n_i = Ae^{-\beta\epsilon_i},A = e^{-\lambda}$
$n_i = N\frac{e^{-\beta\epsilon_i}}{\sum_ie^{-\beta\epsilon_i}}$
$p_i = \frac{e^{-\beta\epsilon_i}}{\sum_ie^{-\beta\epsilon_i}}$
$Z = \sum_{i=1}^{M}e^{-\beta\epsilon_i}$
$T=0,\beta = \infty,p_1=1,p_{i\geq 2}=0$
$T=\infty,\beta=0,p_i = \frac{1}{M}$
$E = \sum_i \epsilon_in_i=N\sum_i\epsilon_ip_i=\frac{N}{Z}\sum_i\epsilon_ie^{-\beta\epsilon_i}=-N\frac{\partial\ln Z}{\partial\beta}$

Constant total energy of system depends on temperature

### Connection

$dQ=TdS$
$$E=U=\sum_i\epsilon_in_i$$
$$dU = \sum_i\epsilon_idn_i+\sum_in_id\epsilon_i$$
因为$\ln \Omega$里面有$dn_i$，所以第一项指的是$dQ$，第二项是$dW$

### Ensembles

$$S_{AB}=S_A+S_B$$

<font color='red'>看看下面三个定义怎么推的</font>

#### Temperature definition

$$\frac{1}{T}\equiv \frac{\partial S}{\partial E}$$

#### Pressure definition

$$\frac{P}{T}=(\frac{\partial S}{\partial V})_E$$

#### Chemical potential definition

$$\frac{\mu}{T}=-(\frac{\partial S}{\partial N})_{E,V}$$

#### Free energy

<font color='red'>看ppt还有tutorial吧</font>

### Partition Function

#### Degeneracy

$Z=e^{-\beta\epsilon_1}+2e^{-\beta\epsilon_2}$

In **equilibrium**, E,S and F are connected t Z by simple equations

$S=\frac{E}{T}+Nk\ln Z$
$F =E-TS$
$S=Nk(\ln Z +T\frac{\ln Z}{\partial T})$

#### Partition function of the system

$Z_N=Z^N$
$Z_N=\sum_ie^{-\beta E_i}$

#### Boltzmann Distribution derivation by system

$p_i=\frac{1}{Z_N}e^{-\beta\epsilon_i}$

## Mixing-Osmosis

### Mixing Entropy and dilute solutions

Lattice models: $N=N_W+N_S$

#### Entropy of Mixing

### Polymer solutions-Mixing entropy

#### Polymer solutions

Entropy of mixing: $\Delta S_{mix}=-k(A\ln \frac{A}{A+B}+B\ln \frac{B}{A+B})$

<font color='red'>看ppt还有tutorial吧</font>

### Osmosis

Entropy of mixing: $\Delta S_{mix}=-k(N_S\ln \frac{N_S}{N_S+N_W}+N_W\ln \frac{N_W}{N_S+N_W})\approx -k(N_S\ln \frac{N_S}{N_W}-N_S)$

Solvation energy = contribution per solute molecule multiplied by the total number of such molecules.
$$dU=\epsilon_s dn_s$$
**First law of thermodynamics**
$\epsilon_sdN_s+\mu_W^0dN_W=TdS+\mu_SdN_S+\mu_WdN_W$
分别令$dN_S,dN_W=0$, 有$\mu_W=\dots, \mu_S = \dots$
使用$\Delta S_{mix}$得到所需偏导
最后得到：$\mu_S=\epsilon_S+kT\ln \frac{n_S}{n_W}$
$$\mu_S=\mu_0+kT\ln \frac{c}{c_0},\enspace c =\frac{n_S}{V},\enspace c_0=\frac{n_W}{V}$$

#### Osmotic Pressure

$$\Delta P=cRT$$
$R = k_BN_A$
if equation contains kB, quantities relate to one particle, otherwise relate to one mole

$$\Delta P=cRT,\enspace c(mol/L)$$

$$\Delta P=ck_BT,\enspace c(个/m^3)$$

### Revere Osmosis

<font color='red'>看tutorial吧</font>

## 5. Binding

### 5.1 Equilibria Binding

#### Ligand-Receptor Binding

<font color='red'>要看一下ppt和tutorial，推$p_{bound}$</font>
Hill function

#### Reactions

$$K_{eq} \equiv \prod_i c_{i-eq}^{\upsilon_i}$$

$\Delta G = \sum_i\upsilon_i\mu_i$
$\mu = \mu_0+kT\ln \frac{c}{c_0}$
$\Delta G =0$
$$\Delta G_0=-k_BT\ln K_{eq}$$

$$\Delta G=k_BT\ln \frac{K}{K_{eq}}$$

### 5.2 Protein folding

<font color='red'>全是叙述，看一下ppt，不看也没问题</font>

## 6. Grand-canonical-cooperativity

### 6.1 Gibbs-Ensemble

finding a given state of the system(charaterized by an energy $E_A$ and nuber of particles $N_A$ is **proportional to the number of states available to the reservoir $R$**) when the system is in this state

$$p_A = \frac{Z}{e^{-\beta(E_A-\mu N_A)}}$$
$$Z=\sum_ie^{-\beta (E_i-\mu N_i)}$$

<font color='red'>看一下tutorial，怎么用公式</font>

### 6.2 Cooperativity

<font color='red'>使用上节公式，看一下tutorial，怎么用公式</font>

## 7. Aminoacids-pKa-ATP-Phosphorylation

### 7.1 Aminoacids-pKa

Proteins, chirality, diversity, peptide bond formation, basic protein structure, from primary to quaternary structure, disease
$K_W=[OH^-][H_3O^+],pH$
$pK_a=pH-\lg \frac{[A^-]}{[AH]}$
得：$\frac{[A^-]}{AH}=\dots$
$$p(A^-)=\frac{1}{1+10^(pK_a-pH)}$$

$$p(HA)=\frac{1}{1+10^(pH-pK_a)}$$

### 7.2 ATP

<font color='red'>看一下tutorial ATP怎么算</font>

### 7.3 Phosphorylation

<font color='red'>看一下tutorial ATP怎么算</font>

## 8. Diffusion

### 8.1 Mean Variance SEM

mean: $\overline{x}=\frac{1}{N}\sum_{i=1}^{N}x_i=\sum_ip_ix_i=\int_{-\infty}^{\infty}xp(x)dx$
$\sigma:$ 三种形式

<font color='red'>推一下ppt概率公式？</font>

### 8.2 Diffusion

Divergence theorem: $\int_V\nabla \cdot \mathbf EdV =\oiiint \mathbf E\cdot d\mathbf A$
Continuity equation: $\nabla\cdot\mathbf j(\mathbf r,t)=-\frac{\partial \rho(\mathbf r,t)}{\partial t}$
Flux: $\mathbf j \equiv \frac{Substance}{Area \times Time} = \rho \mathbf v$
Fick's First Law:
$$j = -D\frac{d\rho}{dx},\enspace \mathbf j=-D\nabla \rho$$

Fick's Second Law:
$$\frac{\partial\rho}{\partial t}=D\nabla^2\rho$$
Solution:
$$\rho(x,t)=\frac{\rho_0}{\sqrt{4\pi Dt}}e^{-\frac{(x-x_0)^2}{4Dt}}$$

$$\rho(\mathbf r,t)=\frac{\rho_0}{\sqrt{4\pi Dt}}e^{-\frac{(\mathbf r-\mathbf r_0)^2}{4Dt}}$$
mean: $\overline{x}(t)=0$
variance: $\sigma_x^2=2Dt$

#### Brownain motion microscopic model: Random walk

$\overline{X}=0$
But what counts is not the mean position, but the spread
$\sigma_x^2=N(\overline{x^2}-\overline{x}^2)=N\overline{x^2}=Nl^2=N\overline{v}l\tau=\overline{v}lt$
又：$\sigma _x^2(t)=2Dt$
$D=\frac{1}{2}\overline{v}l$
Different equations for different situations are in ppt

$$x = \sigma_x = \sqrt{2D}\sqrt{t}$$

#### Diffusion with external forces

$v_D=\mu F=\frac{1}{\gamma}mg$
Down-flux: $j_1=\mu mg\rho$
Up-flux: $j_2=-D\frac{d\rho}{dh}$
$j_1=j_2$
$$\rho(h)=\rho_0e^(-\frac{\mu mg}{D}h)$$

Compared to Boltzman distribution:
$$\rho(h)=\rho_0e^{-\beta mgh},\enspace E=mgh$$

Einstein Relation:
$$D=\mu kT$$

#### Special case for Stoke's Law

Frictional force (drag force) exerted on spherical objects with very small Reynolds numbers (e.g., very small particles)in a continuous viscous fluid
Stoke's Law:
$$F=6\pi R\eta v,\enspace \mu=\frac{v}{F}$$

$$D = \frac{kT}{6\pi R\eta}$$

Flux with external force: $\mathbf{j}=-D\nabla \rho +\mu F\rho$
Diffusion equation with external force: $\frac{\partial \rho}{\partial T}=D\nabla^2\rho-\mu F\nabla\rho$

## 9. Electro

### 9.1 Poission-Boltzmann

#### Electronics for Salty Solutions

Magnetic fields play (almost) no role
$$\nabla \cdot\mathbf{E}(\mathbf{r})=\frac{\rho(\mathbf{r})}{\varepsilon_0}$$
Poisson's equation: $\nabla^2\Phi(\mathbf{r})=-\frac{\rho(\mathbf{r})}{\varepsilon_0}$
Laplace's equation: $\nabla^2 \Phi(\mathbf{r})=0$

$p_{Ion}(\mathbf{r})\propto e^{-\beta E(\mathbf{r})}, \enspace E(\mathbf{r})=q\Phi(\mathbf{r})$
**Charge density of the ion cloud**:
$$\rho_i(\mathbf{r})=q_in_i^{\infty}e^{-\beta q_i\Phi(\mathbf{r})}$$

**total charge density**:
$$\rho_{Ions}(\mathbf{r})=\sum_iq_in_i^{\infty}e^{-\beta q_i\Phi(\mathbf{r})}$$

**Poisson-Boltzman equation**
$$\nabla^2\Phi(\mathbf{r})=-\frac{e}{\varepsilon}\sum_{i=1}^nz_ic_{i,\infty}e^{-\beta z_ie\Phi(\mathbf{r})}-\frac{\rho(\mathbf{r})}{\varepsilon}$$

<font color='red'>看一下tutorial</font>

### 9.2 Debye-Huckel theory

#### Bjerrum length

Balance the ekectrostatic interaction energy:
$$\frac{e^2}{4\pi \varepsilon_0Dl_B}=k_BT$$

$$l_B=\frac{e^2}{4\pi\varepsilon_0Dk_BT}$$

#### Debye-Huckel theory

**Debye-Huckle equation:**
$$\nabla^2\Phi(\mathbf{r})=\frac{e^2\beta}{\varepsilon}(\sum_iz_i^2n_i^{\infty})\Phi(\mathbf{r})$$

**Debye length**
$$\lambda_D=(\frac{\varepsilon k_BT}{e^2\sum_iz_i^2n_i^{\infty}})^{\frac{1}{2}}$$

**Ionic strength**
$$I\equiv \frac{1}{2}\sum_iz_i^2n_i^{\infty}$$

$$\nabla^2\Phi(\mathbf{r})=\frac{1}{\lambda_D^2}\Phi(\mathbf{r})$$

<font color='red'>看一下tutorial中公式怎么用的</font>

## 10. Polymers-CentralLiit

### 10.1 Intro-Polymers

### 10.2 Binomial distribution

$$
p(k) = C_n^k p^k(1-p)^{n-k}
$$

Special case p=0.5

Expectation value:

If n=1, p(1)=p.
$\langle x \rangle = \sum_{i=1}^1 1\cdot p = p$
$\langle x+y \rangle = \langle x \rangle \langle y \rangle $ 
therefore $\langle x \rangle = np$

#### Poisson limit theorem

n较大, p较小时，二项分布可以近似为泊松分布

### 10.3 Polymers-Random walk Central limit theorem

Atoms don't always matter.

#### Random walks: freely jointed chain (FJC)

几段等长的segment连接，每一段都可以是任意方向，each segment 的长度叫做 **Kuhn length (a)**

##### Chain: 1D random walks

The expected value of the walker's distance from the origin, R, after N steps is 0.
**Similar to diffusion**
mean = 0, variance = $Na^2$
**R为end-to-end distance**

#### Central limit theorem

For a set of N random variables {x_n} with finite mean and variance, the sum X of {x_n} will tend towards a **Gaussian** distribution *regardless of the distribution of x_n*

### 10.4 Polymers - Random walk

Sharpness of the Gause curve:
$$
\sigma_x = \sqrt{N}
$$

$$
\frac{\sigma_x}{N} = \frac{1}{\sqrt{N}}
$$

As N becomes larger, the RMS distance increases, but the relative (to N) RMS distance becomes tiny, so it is extremely unlikely to be far (as compared to N) from the pub door.

#### Polymer end‐to‐end vector (1D)

$\langle x^2 \rangle = Na^2$
$p(x,N) \frac{1}{2\pi Na^2}e^{-\frac{x^2}{2Na^2}}$

#### 3D case

$$
P(R;N) = P(R;N)(1D case)^3
$$

Similar to diffusion

##### End-to-end distance R

$P_{3d}(N,R)4\pi R^2dR$
**3D Maximum of distribution NOT at zero! Different from the 1D case.**

##### Limitation of Guassian model

Wrong for $R>Bb$

$P(R;N)$ is sharply peaked at R=0

## 11. FJC-Kuhn_Rg

### 11.1 FJC

$C_n=\frac{1}{n}\sum_{i=1}^nC_i'$: Flory's characteristic ratio
$\langle R^2 \rangle = l^2\sum_{i=1}^nC_i'=C_nnl^2$

$C_{\infty}$: Flory's characteristic ratio, for large chain

#### Freely rotating chain model (FRC)

different from FJC: $\theta$ is const. same for all

$$
\langle R^2 \rangle = nl^2 + 2nl^2\frac{\cos \theta}{1- \cos \theta} = nl^2\frac{1+\cos \theta}{1-\cos \theta}
$$

$$
C_{\infty} = \frac{1+\cos \theta}{1-\cos \theta}
$$

#### Persistence length

$$
s_p = -\frac{1}{\ln(\cos \theta)}
$$

**For rigid model:** $L_p=a$
**For Continuous model:** tangent correlation function
$$
\langle t(s)\cdot t(u) \rangle = e^{-|s-u|/\zeta _p}
$$
**For Worm-like chain(WLC) modle：DNA:**
small $\theta, \ln(\cos \theta)\approx -\frac{\theta ^2}{2}$
$$
l_p = s_pl = l\frac{2}{\theta ^2}
$$

$$
C_{\infty} = \frac{1+\cos \theta}{1-\cos \theta} \approx \frac{4}{\theta ^2}
$$

$$
b=l\frac{C_{\infty}}{\cos (\frac{\theta}{2})}\approx l \frac{4}{\theta ^2}=2l_p
$$

Kuhn lebgth b is twice the persistence length.

### 11.2 Polymers - Radius of Gyration(回转半径)

**The radius of gyration is used because it can be easily measured experimentally.**

Radisu of gyration determined via scattering experiments

**Tutorial 证明**
<font color='red'>看一下tutorial中公式怎么用的</font>

## 12. EntropicForce-AFM

### 12.1 Entropic Force

$$
\Omega (N,x) \propto p(N,x), S = k_B\ln \Omega
$$

$$
\Delta S = -\frac{k_B}{2Na^2}x^2
$$

$$
\Delta G = \Delta U - T\Delta S
$$

The ideal chain has no internal energy U

$$
\therefore \Delta G = -\frac{k_BT}{2Na^2}x^2
$$

$$
f = \frac{\partial G(N,x)}{\partial x} = \frac{k_BT}{Na^2}x
$$

Hooke's law:
entropic spring constant
$$
k = \frac{k_BT}{Na^2}
$$

It's easier to stretch polymers with:

- large numbers of monomers N
- large monomer size a
- at lower temperature T

*The entropic nature of elasticity in polymers distinguishes them from other materials.*

**The ideal chain can be thought of as an entropic spring and obeys Hooke's law for elongations much smaller than the maximum elogation.**
可以把上面那个x换成end-to-end vector R

<font color='red'>等会看ppt吧，写不下去了</font>

### 12.2 DNA Properties and genetics

structure and packing
**Polymerase chain reaction**

## 13. DNA-Genetics-Bending

### 13.1 DNA-Genetics

GENETICS NETWORKS:
DOING THE RIGHT THING AT RIGHT TIME
<font color='red'>maybe 再看一下ppt上概率那一张</font>

### 13.1 Bending

Three modes of deformation of a beam:

- stretching
- bending
- twisting

if a bond stretches by $\Delta a$, the beam stretches by $\Delta L = \frac{L_0}{a_0}\Delta a$

$$
\epsilon = \frac{\Delta L}{L_0} = \frac{\Delta a}{a_0}
$$

<font color='red'>看ppt上公式吧，还有tutorial</font>

## 14. DNA-Looping-packing-Viruses

### 14.1 DNA-Looping

$$
E_{bend} = \frac{EIL}{2R^2}, L = 2\pi R
$$

$$
\therefore E_{loop} = \frac{\pi EI}{R} = \frac{2\pi^2EI}{L}, \zeta_p / x_p = \frac{EI}{k_BT}
$$

$$
E_{loop} = k_BT2\pi^2\frac{x_p^2}{L} = k_BT2\pi^2\frac{x_p^2}{d\cdot N_{bp}}, L = d \cdot N_{bp}
$$

$$
\Delta S_{loop} = k_B \ln p_o = k_B(-\frac{3}{2}\ln N_{bp}+const)
$$

$p_o$ is the probability of loop formation, $p_o \propto N_{bp}^{-\frac{3}{2}}$

$$
\Delta G_{loop} = \Delta E_{loop} - T\Delta S_{loop}
$$

### 14.2 DNA packing, Viruses

Viruses as Charged Spheres
$q = -2e$ per base pair

$$
U_{el} = \int dU = \int _0^R V(r)dq
$$

**DNA in solution:** $\lambda _D$
**DNA in capsid**

$$
W = NkT \ln (\frac{V_{cloud}}{V_{capsid}})
$$

<font color='red'>ppt上单分子技术</font>

## 15. Reactions

<font color='red'>看一下tutorial中公式怎么用的或者pdf推导</font>

Some examples

***Cytoskeleton polymerization***

$$
\frac{dn}{dt} = k_{on}(c_)-\frac{Mn(t)}{V}-k_{off}
$$

## 16. Membrane-Surfaces

### 16.1 Membranes-Intro

- 2D fluid: Lateral diffusion
- MP can also diffuse
- Lipid flip-flop

Lipids self assemble

The effective shape of a lipid molucule is described by a packing parameter:

$$
P = \frac{v}{al}
$$

Membrane Permeability

$$
j = P\Delta c
$$

P为Permeability coefficient

All membranes indergo spontaneous shape changes and fluctuations due to thermal energy or application of external forces

Membrane fusion and budding

### 16.2 Surfaces Math Background

<font color='red'>看一下tutorial中公式怎么用的, ppt 要是有时间也可以看看吧</font>

$$
dE_{surface} = [\frac{\kappa}{2}(\frac{1}{R_1}+\frac{1}{R_2}-\frac{2}{R_0})^2+\frac{\kappa G}{R_1R_2}]dA
$$

### 16.3 Membrane deformation and curvature

**bending**
curvature: height function

$$
\kappa = \frac{1}{R} = \frac{\partial ^2h}{\partial x^2}
$$

对角化Hessian matrix, 本征值即为两个方向上的曲率

$$
G_{bend} = \frac{K_b}{2}\int da [\kappa _1+\kappa _2]^2
$$

$$
G_{thickness} = \frac{K_t}{2} \int (\frac{w-w_0}{w_0})^2da
$$

$$
G_{stretch} = \frac{K_a}{2} \int (\frac{\Delta a}{a_0})^2 da
$$

## 17. Membrane-Potential

### 17.1 MPs

#### Membrane proteins

#### MP structure determination

### 17.2 Membrane Potential

$$
c(x) \sim e^{\beta E(x)} = e^{-\beta z eV(x)}
$$

$$
\Delta V = V_2 -V_1 = \frac{k_BT}{ze} \ln \frac{c_1}{c_2}
$$

$$
\Delta \epsilon = \Delta \epsilon _{conf} - p \frac{V_{mem}}{d}
$$

$$
p_{open} = \frac{e^{-\beta \Delta \epsilon}}{1+e^{-\beta \Delta \epsilon}}
$$

Nerst equation relates chemical potential to electric potential

## 18. Membrane-PatchClamp

### 18.1 Vesicles

For a vesicle of radius R:

$$
GG_{bend} = 8\pi K_b
$$

For a Cylinder:

$$
G_{bend} = K_b\pi \frac{L}{R}
$$

### 18.2 Surface Tension

表面处：$E = -\frac{z}{2}\epsilon$
内部：$E = -z\epsilon$

$$
\Delta E = +\frac{z}{2}\epsilon
$$

Surface tension $\gamma$: Energy per unit area of the surface = surface energy density

$$
\gamma = \frac{F}{2L} = \frac{dW}{dA}
$$

$dA$ is total area change, 2 sides

**Laplace-Young law:**

$$
\Delta P = \frac{2\gamma}{R}
$$

<font color='red'>看tutorial怎么用</font>

### 18.3 Membranes: Patch clamp

<font color='red'>看tutorial怎么用</font>

