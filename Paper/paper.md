---
title: '`INS_IEP`: A MATLAB package for fitting peaks of Inelastic Neutron Scattering data in spin cluster systems'
tags:
  - MATLAB 
  - Inelastic Neutron Scattering
  - Inverse Eigenvalue Problems
  - Deflation
  - NonLinear Least Squares
authors:
  - name: Alban Bloor Riley
    orcid: 0000-0001-9431-709X
    affiliation: 1
    corresponding: true
  - name: Marcus Webb
    orcid: 0000-0002-9440-5361
    affiliation: 1
  - name: Michael L. Baker
    orcid: 0000-0002-8246-3177
    affiliation: 2
affiliations:
  - name: The University of Manchester, Department of Mathematics
    index: 1
  - name: The University of Manchester, Department of Chemistry
    index: 2
date: June 2025
bibliography: paper.bib

---

# Summary

Inelastic neutron scattering (INS) is a spectroscopic technique that can be used to measure the magnetic excitations in materials with interacting electron spins. For samples composed of finite-size clusters of magnetic moment-carrying atoms, such as single ions or molecular-based magnets, INS experiments yield inelastic excitation with energies that correspond to the eigenvalues of the spin Hamiltonian of the material being studied. Fitting a model spin system to these experimental eigenvalues can be formulated as an inverse eigenvalue problem (IEP), where the matrix formed is the spin Hamiltonian operator of the sample molecule. Solving this IEP is less computationally expensive than fitting the full INS data, so it can be used as an initial proxy for the full fitting problem. `INS_IEP` is a MATLAB package that uses deflated numerical optimisation methods to systematically find multiple solutions to this IEP. The package requires and is fully compatible with `easyspin` [@stoll_easyspin_2006], a package for solving fitting problems in electron paramagnetic resonance (EPR).

# Statement of need

Neutrons are an excellent bulk probe of material properties since they carry no charge and therefore penetrate deeply into matter. Neutrons also carry a quantum spin of a half, making them a sensitive probe of magnetism [@squires_introduction_2012]. Reactors and spallation sources with dedicated high-flux neutron sources serve the international research community with neutron scattering experiment capabilities for material research.
Inelastic neutron scattering (INS) is one such experimental technique that can be used to study magnetism. In an INS experiment, a sample under investigation is irradiated with a beam of neutrons and the scattered neutron energy and momentum transfer are detected. For samples composed of finite-size clusters of magnetic moment-carrying atoms, such as single ions or molecular-based magnets, the detected neutron energy transfer gives direct access to the quantum spin excitations [@furrer_magnetic_2013;@baker_neutron_2012;@baker_spectroscopy_2014]. 

The energy of such excitations relates to the energy difference between eigenvalues of the Hamiltonian matrix that describes the quantum spin dynamics of the compound in question. Single-ion and molecular-based magnets are studied as prototype components (quantum bits, sensors) for quantum technologies. INS can therefore provide crucial information concerning the precise quantum properties of such systems. However, to relate the INS experimental results to the Hamiltonian that describes quantum spin dynamics requires parameterisation of matrix elements such that a set of eigenvalues and eigenstates matching the experiment is determined. This situation is known as the inverse eigenvalue problem. 

To date, this problem is addressed in an iterative process where parameters of the Hamiltonian are varied manually, often one at a time, and the resultant eigenvalues compared to the experimental values - each such iteration requires an eigendecomposition of the Hamiltonian matrix. INS_IEP presents an elegant solution to solving this problem, using algorithms to calculate multiple parameter sets that minimise the difference in eigenvalues, reducing the number of Hamiltonian matrix diagonalisations, and providing a more robust method to reliably extract an accurate spin Hamiltonian model from INS experimental data. 

## Current State of the Art 

There currently exists fitting software designed to solve this problem, such as:

- `SPECTRE` [@a_t_boothroyd_spectre_1990] uses the NAG Fortran Library for minimisation and matrix diagonalisation.
- `FOCUS` [@fabi_focus_1995] uses Monte-Carlo and general quasi-Newton methods to find the parameters.
- `PyCrystalField` [@scheie_pycrystalfield_2021] -  fitting uses the SciPy minimize library. 

All of these methods, however, only work for single ions (i.e., no exchange coupling) and only calculate one minimising system. On the other hand, `INS_IEP` can fit INS data from samples with multiple spin centres, and it can also find multiple minimising systems without multiple initial guesses to the parameter set.

# Key Concepts

## The Spin Hamiltonian

The Spin Hamiltonian, $H$, is an approximation of the Hamiltonian that uses spin coordinates instead of orbital coordinates, and is widely used to model data arising from many spectroscopy techniques [@launay_electrons_2014]. It can be modeled as a linear combination of interaction terms; in this package, we will use the  zero field interaction, $H_{ZFI}$, and the electron-electron interaction, $H_{EEI}$:

$$H = H_{ZFI} + H_{EEI}.$$

Both of these terms can themselves be modelled as the linear sum of other basis matrices. The zero-field interaction can be written as:

$$H_{ZFI} = \sum_{-k\leq q \leq k} B^q_kO^q_k$$

where the $O^q_k$ are Stevens Operators [@rudowicz_generalization_2004], and $B^q_k$ the associated parameter. When there are multiple spin centres, it is necessary to take Kronecker products of the operator with identity matrices of the appropriate size for each spin centre.

When there are multiple spin centres, it is also necessary to include an electron-electron interaction term, $H_{EEI}$. This term will be the sum of interaction terms between each pair of spin centres:

$$H_{EEI} = -\sum_{i\neq j} J_{ij}  S_i\cdot  S_j$$

where $S_i$ is the vector of spin operators $S_i = [S_x, S_y, S_z]$ for the $i$-th spin centre, and $J_{ij}$ is the parameter to be found that represents the strength of interaction between the two spin centres. Note that in the isotropic case, $J$ can be thought of as a scalar value, but in the anisotropic case will be a matrix where the off-diagonals are skew-symmetric (often the off-diagonals are assumed to be zero). While the summation is in theory over all spin centre combinations, in practice many of these contributions will be negligible - often only the nearest neighbour interactions are significant. 

It is important to mention that these matrix operators can be very large. The size is defined by the number of spin centres ($n$) and the spin $(S_i)$ of each spin centre. The dimension of matrices is given by:
$$ \prod_i^n (2S_i+1).$$
However, the operators are highly sparse, which means that it is possible to use eigensolvers that can take advantage of this sparsity.


## Inverse Eigenvalue Problem


The INS experiments provide eigenvalues of the Spin Hamiltonian matrix of the sample. The task of calculating the matrix from the eigenvalues is an inverse eigenvalue problem:

Let $A(x)$ be the affine family of matrices,

$$A( x) = A_0 + \sum^\ell_{i=1} x_i A_i,$$

where $x\in\mathbb R^\ell$ and $A_0,\dots,A_\ell \in \mathbb C^{n\times n}$ are linearly independent Hermitian matrices, and denote the ordered eigenvalues of $A(x)$ as $\lambda_1(x)\leq\dots\leq\lambda_n(x)$.
Then the least squares inverse eigenvalue problem (LSIEP) is to find the parameters $x \in \mathbb R^\ell$ that minimise:

$$F(x) = \frac 1 2 ||r(x)||^2_2 = \frac 1 2 \sum^m_{i=1}(\lambda_i(x) - \lambda_i^*)^2$$

where $\lambda^*_1\leq\ldots\leq\lambda_m^*$ are the experimental eigenvalues [@chu_inverse_2005]. In the case of INS, fitting the $A_i$ basis matrices will be a combination of Stevens operators and electron-electron exchange terms. The IEP described above is formulated as a least squares problem because the number of eigenvalues that can be probed by INS experiments is often a small subset of the full spectrum. Due to the low temperatures at which these experiments are performed (which can be as low as 1K), it is generally the smallest eigenvalues that are involved. Note also that since it is the energy difference between the eigenvalues that is probed we actually have to modify the IEP - either by adding an additional parameter (an identity matrix) that shifts the values of the eigenvalues, or by changing the above formula for $F$ to directly sum the difference in eigenvalues thereby reducing the number of residual equations in $r(x)$ by one. 

As far as we are aware, this is the first time that the fitting of INS data has been explicitly formulated as an IEP. An advantage of this formulation is that there are explicit formulas for the derivatives of $r(x)$. The first derivative (Jacobian) is:

$$ J_r(x) = \begin{pmatrix}
        q_1(x)^TA_1q_1(x)&\dots &q_1(x)^TA_\ell q_1(x)\\
        \vdots&\ddots&\vdots\\
        q_m(x)^TA_1q_m(x)&\dots& q_m(x)^TA_\ell q_m(x)
    \end{pmatrix},$$    

and the second derivative (Hessian) is:

$$ (H_{r})_{ij}   = 2\sum^m_{k=1} (\lambda_k-\lambda_k^*) \sum^m_{\substack            {t=1\\\lambda_t\neq\lambda_k}} \frac{(q_t^TA_iq_k)(q_t^TA_jq_k)}{\lambda_k-\lambda_t}.
$$
Another advantage is that the number of constraints to fit is much smaller than fitting the spectrum itself, as it corresponds to fitting only the locations of the peaks of the spectrum.

## Methods

All of the methods used are iterative schemes of the form $x^{k+1} = x^k +p^k$ where the step $p^k$ uniquely defines each algorithm:

-  Newton's method: $p^k = (J_r^TJ_r + H_rr)^{-1}J_r^Tr$ [@nocedal_numerical_2006]
-  Gauss-Newton method: $p^k = (J_r^TJ_r)^{-1}J_r^Tr$ [@nocedal_numerical_2006]
-  Lift and Projection Method: $p^k = B^{-1}J_r^Tr$ [@bloor_riley_riemannian_2025]

Where the matrix $B$ is the Gram matrix formed from the Frobenius inner products of the basis matrices: $B_{ij} = \langle A_i, A_j\rangle_F$. The Lift and Projection method is a Riemannian Gradient descent method [@bloor_riley_riemannian_2025], inspired by the Lift and Projection method [@chu_inverse_2005], specifically designed for solving IEPs. In [@bloor_riley_riemannian_2025], it is proven that the method is a strictly descending algorithm and reduces the value of the objective function every step. Both the deflated Gauss-Newton method and the Riemannian Gradient descent Lift and Projection method are new methods designed for this package [@bloor_riley_deflation_2025;@bloor_riley_riemannian_2025]. 

The $m$ eigenvalues required to calculate $J_r$ are calculated using MATLAB's `eigs` command, which is an efficient method based on Krylov subspaces to calculate a small subset of the eigenvalues of the matrix. This is invaluable in this setting because computing the full set of eigenvalues by exact diagonalisation is, in some cases, infeasible [@bloor_riley_riemannian_2025].

## Deflation 

The number of eigenvalues that can be probed via INS experiments varies depending on the equipment and sample in question, meaning that the fitting problem is often under- or even over-determined. The IEP is also highly nonlinear and, due to the experimental nature of the data, may be ill-posed. One consequence of this is that the solution space may be very 'bumpy', and there may exist many local minimisers to the problem. For example, in \autoref{fig:contour}, there are clearly 4 distinct solutions (for more details, see Example 1 and the file Example1_Mn12.m in the examples folder). We seek to solve the problem of multiple local minima by the use of Deflation, a numerical technique used to find multiple solutions to systems of equations [@farrell_deflation_2015]. Fortunately, it is cheap to apply deflation for the above methods, and it is simply a change to the length of the step - notably, this means that the direction of each step does not change. It is proven in [@bloor_riley_deflation_2025] that the deflated methods will not converge to deflated points.  The usual requirements still apply to the convergence of the new methods - that the initial guess is close enough to the new minimum, and that the Jacobian is full rank in a neighbourhood around that minimum. The rate of convergence of the deflated methods is also more complicated, although the number of iterations required to converge can go up with the number of deflations; this is not a strict correlation, as can be seen in  \autoref{fig:convergence}.

![Contour plot of how $F$ varies with the two parameters $B^2_2$ and $B_4^4$ for the molecule Mn\_12 as described in Example 1. There are four locally minimising parameter pairs corresponding to the four blue regions \label{fig:contour}](Mn12_contour.png){ width=90%}

![Comparison of the convergence behaviour for computing each solution in Example 1. The gradient of the lines as the method approaches the solution shows how the methods are quadratically convergent local to a minimum.  \label{fig:convergence}](Mn12_convergence.png){ width=90%}

# Examples
## Example 1 - Mn12 
The first example we will look at is Manganese-12-acetate. This is a well-known example in the INS and magnetism community, as one of the first molecules that behaves like a nano-sized magnet with a molecular magnetic coercivity as well as its role in the research of quantum tunnelling of magnetisation [@friedman_macroscopic_1996;@sessoli_magnetic_1993]. 

The Spin Hamiltonian of this system, using the giant spin approximation, can be represented as a $21\times21$ matrix modelled using 4 Stevens operators [@bircher_transverse_2004]:

$$H = B^0_2O^0_2 + B^0_4O^0_4 +B^2_2O^2_2 +B^4_4O^4_4 \in \mathbb R^{21\times 21} $$

We utilise the same spin system syntax as `easyspin`, so to set up the problem, we first set up the model, along with initial guesses for the parameters:
```matlab
  %One spin centre (because giant spin approximation)
  Sys0.S=10; 
  %Four Stevens operators 
  Sys0.B2 = [-100,0,-1000,0,0];   
  Sys0.B4 = [-1,0,0,0,-1,0,0,0,0];
```
Then we input the experimental eigenvalues - these are typically shifted such that the smallest eigenvalue is zero - and define which parameters to fit. Note that all values given must be in gigahertz, so it may be useful to use conversions.
```matlab
  rcm = 29979.2458;   meV = rcm*8.065;  %Conversions values
  %Input calculated eigenvalues:
  Exp.ev = [0,0,1.24,1.24,2.3,2.3,3.18,3.18,3.91,3.91,4.5,4.5,
           4.97,4.97,5.32,5.32,5.54,5.59,5.69,5.75,5.78].*meV;
  %Note that these eigenvalues are simulated from the parameters given in [@bircher_transverse_2004].

  %Vary all non-zero parameters (no Fixed parameters):
  Vary = Sys0; 
```
Then all that is required is to call `INS_IEP` with these three inputs:
```matlab
  SysOut = INS_IEP(Sys0,Vary,Exp);
```
If we wish to find all four solutions as shown in \autoref{fig:contour}, then we use the additional option:

```matlab
  Opt.NDeflations = 4;
  SysOut = INS_IEP(Sys0,Vary,Exp,Opt);
```
In this case, SysOut will be an array of four spin structures, each containing a distinct locally optimal solution. It is possible to access information about the convergence of each deflation by using `SysOut.Output`. For example, by utilising the iterates recorded, stored in `SysOut.Output.Iterates`, it is possible to plot a graph of convergence, as can be seen in \autoref{fig:convergence}. The Output structure also contains the value of $F$ at the final point, as well as the number of iterations it took to get there.

A full list of options is provided in the help of INS_IEP.

## Example 2 - Chromium(iii) Horseshoes
The second example concerns an antiferro-magnetically coupled chain of six chromium(III) ions [@baker_varying_2011]. Because there are multiple spin centres, an electron-electron interaction term is required. The spin Hamiltonian is a 4096×4096 matrix composed of two Stevens operators and one interaction term, since it is known a priori that each spin centre will have the same value parameters. We pin the parameters here by setting the initial guess to the same value:

```matlab
  Sys0.S = [1.5 1.5 1.5 1.5 1.5 1.5];
  Sys0.B2 =  [1     0    -1     0     0;
             1     0    -1     0     0;
             1     0    -1     0     0;
             1     0    -1     0     0;
             1     0    -1     0     0;
             1     0    -1     0     0];
  Sys0.J = [100,0,0,0,0,100,0,0,0,100,0,0,100,0,100];
  Vary = Sys0;
  Exp.ev = [0,0.355,0.457,0.497,1.576,1.577,1.592,1.629,1.632,
            2.97,2.98,3.002,3.004,3.01,3.038,3.821,3.824,3.827,
            3.837,3.856,3.879,3.888,3.895,3.903].*meV;
```
Note that only 24 eigenvalues were found experimentally, so this will form a partial LSIEP. To find the solution system is as simple as:

```matlab
  Opt.GradientTolerance = 1e-3; %Use additional stopping criterion.
  SysOut= INS_IEP(Sys1,Vary1,Exp,Opt);
```
It is possible to find multiple minimising systems even if they do not make any sense physically. However, due to the scaling of the problem, a change in the default deflation parameters is necessary:
```matlab
  Opt = struct('NDeflations',5,'Sigma',1e-7,'StepTolerance',1e-3 );
  SysOut= INS_IEP(Sys0,Vary,Exp,Opt);
```
The output contains five different spin systems that all have the same eigenvalues as input. Like in example one, some of the minima found differ only by a sign change of the $B^2_2$ parameter, but all minima have the same exchange term. 

Additional examples can be found in the Examples folder. 

Note also that using `mint` [@baker_mint_2022], which is fully compatible with `INS_IEP`, it is possible to simulate the INS spectrum of any calculated system, which can then be compared to the experimental spectrum. A description of how to do this can also be found in the examples folder. 

# Acknowledgements

We thank the reviewers, Garrett Granroth and Chen Zhang, for their help in improving the paper. ABR thanks the University of Manchester for a Dean’s Doctoral Scholarship. MW thanks the Polish National Science Centre (SONATA-BIS-9), project no. 2019/34/E/ST1/00390, for the funding that supported some of this research. 

# References

