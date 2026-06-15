# **SMPNP - a size modified Poisson Nernst Planck simulation tool**

<!-- [![DOI](https://zenodo.org/badge/DOI.svg)](https://doi.org/your-doi)  -->
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  

## **Overview**
This repository contains the source code for the research project **Nanoporous Carbon Coating of Separator Boosts Rate Capability of Cathodes in Lithium-Ion Batteries**, which was published in **Advanced Materials**. The code refers to the numerical simulation portion of the manuscript. Specifically, it **solves a generalized Poisson Nernst Planck system in a two-dimensional nanoslit geometry using the VoronoiFVM.jl package.**

\ud83d\udcdd **Paper DOI**: [Insert DOI]  
\ud83d\udcc1 **Data DOI (if available)**: [Insert DOI for dataset]  
---

## **Installation**
### **Prerequisites**
Before running the code, ensure you have the following installed:
- **Julia**
- **Python 3.10.16** (recommended: `pyenv` or `conda`)
- **Git** for cloning the repository
- **Required libraries** (see `requirements.txt`)

### **Setup Instructions**
#### **Using `pip`**
```sh
git clone [repo-url]
cd [repo-name]
python -m venv venv  # Create a virtual environment
source venv/bin/activate  # Activate it (use `venv\Scripts\activate` on Windows)
pip install -r requirements.txt  # Install dependencies
```

#### **Using `conda` (alternative)**
```sh
conda env create -f environment.yml
conda activate [env-name]
```

#### **Using `pyenv`**
If you use `pyenv` for managing Python versions, run:
```sh
bash pyenv_setup.sh
pyenv activate my_project_env
```

#### **Using Docker (optional)**
For reproducibility, you can use the provided Docker container:
```sh
docker build -t my_project .
docker run --rm -it my_project
```
---


## **Governing equations & numerics**
### **Primary variables**
The code solves the coupled Poisson-Nernst-Planck equations for the following quantities:
$$c_- \quad \text{concentration of anions}$$
$$c_+ \quad \text{concentration of cations}$$
$$\phi \quad \text{electrical potential}$$
### **Governing equations**
The governing equations for the generalized (size-modified) Poisson Nernst Planck are:

$$\frac{\partial c_i}{\partial t} = \nabla \cdot (M_i c_i \nabla \mu_i)\\
\nabla^2\phi= -\frac{1}{\varepsilon}\sum_{j = 1}^{2} z_jec_j,
$$

where $i = 1,2 = \text{Li}^+,\text{PF}_6^-$, $c_i$ is the concentration (mol/m$^3$) of the $i^{\text{th}}$ ion, $t$ is time (s), $M_i = D_i/(k_{\text{B}}T)$ is the ion mobility (m$^2$/(Js)) where $D_i$ is the diffusion coefficient (m$^2$/s), $k_B$ is the Boltzmann constant, $T$ is the temperature (K), $\varepsilon$ is the permittivity (F/m), $\phi$ is the electric potential (V) and $\mu_i$ (J) that accounts for the ion size effect is given by 

$$
\mu_i = k_{\text{B}}T \ln(c_ia_i^3) -\frac{k_{\text{B}}T a_i^3}{a_0^3}\ln\left(1-\sum_{j=1}^{2}c_ja_j^3\right) + z_ie\phi,
$$
with $m_i = a_i^3c_{\text{ref}}$ and $m_0 = a_0^3c_{\text{ref}}$, where $a_i$ is the excluded size (m) of the $i^{\text{th}}$ ion and $a_0$ is the excluded size (m) of the solvent molecule.

### **Non-dimensionalization**
The variables in governing equations are dimensionless. They are made non-dimensional by the following scalings wich can also be found in [^bazant_lecture_notes]

$$\tilde{\mathbf{x}} = \frac{\mathbf{x}}{L_{\text{ref}}}, \quad \tilde{c_i} = \frac{c_i}{c_{\text{ref}}}, \quad \tilde{\varepsilon} = \frac{\varepsilon}{\varepsilon_{\text{ref}}},\quad \tilde{D_i} = \frac{D_i}{D_{\text{ref}}}, \quad \tilde{t} =\frac{t}{\left(L_{\text{ref}}^2/D_{\text{ref}}\right)},\quad \tilde{\phi} = \frac{e\phi}{k_BT}, \tilde{\mathbf{\nabla}} = L_{\text{ref}}\mathbf{\nabla}$$

Here, $\mathbf{x} = (x,y)^T$ is the two-dimensional space vector, $L_{\text{ref}}$ (m) is a chosen reference length of the system, $c_{\text{ref}}$ (mol/m$^3$), $D_{\text{ref}}$ (m$^2$/s), $\varepsilon_{\text{ref}}$ (F/m) are reference scales for concentration, diffusion coefficient and permittivity with the same units as the non-reference quantities.

### **Dimensionless governing equations**
The dimensionless governing PNP equations then read

$$\frac{\partial \tilde{c}_i}{\partial \tilde{t}} = \tilde{\nabla} \cdot (\tilde{D}_i \tilde{c}_i \tilde{\nabla} \tilde{\mu}_i)\\
 \tilde{\lambda}_D^2\tilde{\nabla}^2\tilde{\phi}= -\frac{4\pi}{\tilde{\varepsilon}}\sum_j z_j\tilde{c_j}
$$

where $\tilde{\lambda}_D = \lambda_D/L_{\text{ref}}=\sqrt{\varepsilon_{\text{ref}} k_BT/(\sum_i (z_ie)^2 c_{\text{ref}})}/L_{\text{ref}}$ is the Debye length. The relation for the chemical potential $\tilde{\mu}_i$ is given as 

$$
\tilde{\mu}_i = \ln\left(\tilde{c}_i m_i\right) - \frac{m_i}{m_0}\ln\left({1-\sum_j\tilde{c}_jm_j}\right) + z_i\tilde{\phi}.
$$

## Numerics
The dimensionless governing equations, together with the boundary conditions described in the paper, are solved until steady-state using a backward Euler finite volume method. We implemented it on top of the open-source package [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl). The [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl) package uses Voronoi meshes for the discretization in space and a two-point flux approximation for computing fluxes between different control volumes. The non-linear system of equations is solved with a Newton-Raphson algorithm, with the Jacobian being computed using automatic differentiation.

## Code Verificaiton
To verify the implementation, we compare the numerically computed solutions for two cases against analytical and semi-analytical expressions.

### Analytical PB distribution
For the classical Poisson-Nernst-Planck equations where the ion sizes are not taken into account, there is an analytical solution for the double layer profile. We setup a one-dimensional problem from $x = 0$ to $x = L$ with the following boundary conditions

$$\phi = \tilde{\phi}_{\text{wall}}, J_+ = J_-= 0\qquad \text{at } x = 0$$
$$\phi = 0, c_+ = c_-= 1\qquad \text{at } x = 1$$

All other parameter can be found in the [input file](./test/verification_simple_pnp.yml) for that test case .

The analytical solution is given by

$$ \phi(x) = 2\log\left(\frac{a+b\exp{\left(-x/\lambda_D\right)}}{a-b\exp{\left(-x/\lambda_D\right)}}\right), \qquad a = \exp{\left(\frac{\phi_{\text{wall}}}{2}\right)}+ 1,\quad b = \exp{\left(\frac{\phi_{\text{wall}}}{2}\right)}- 1 $$

$$c_\pm(x) = \exp\left(\mp \phi(x)\right)$$



## Concentration distribution in a size-modified case (dirichlet)
For the case with finite ion sizes, we setup a first test case that still uses the dirichlet boundary conditions for potential and concentration at $x = 1$. The chosen boundary conditions for this case are given as 

$$\phi = \tilde{\phi}_{\text{wall}}, J_+ = J_-= 0\qquad \text{at } x = 0$$

$$\phi = 0, c_+ = c_-= 1\qquad \text{at } x = 1$$

which are the same boundary conditions as for the simple PNP case above. The difference here is that we now solve the size-modified PNP equations that take the finite ion sizes into account. These ion sizes are given as

$$a_+ = a_- = a_0 = 5Angstrom$$

All other parameter can be found in the [input file](./test/verification_finite_size_pnp_dirichlet.yml) for that test case.

For this problem, we have a semi-analytical solution. Given the numerical solution for the potential $\phi(x)$ as well as the two dirichlet values for the concentrations at $x = L$, the solutions for the concentrations $c_+(x)$ and $c_-(x)$ are given as follows

$$ c_-(x) = \frac{1}{m_-}\left(\frac{A_- - \frac{A_-A_+}{1+A_+}}{1+A_- - \frac{A_-A_+}{1+A_+}}\right)$$

$$ c_+(x) = \frac{1}{m_+}\left(\frac{A_+ - \frac{A_-A_+}{1+A_-}}{1+A_+ - \frac{A_-A_+}{1+A_-}}\right)$$

$$A_\pm = \exp\left(\left[\Phi_\pm \mp \phi(x)\right]\frac{m_\pm}{m_0}\right)$$

$$\Phi_\pm = \frac{m_\pm}{m_0}\log\left(\frac{\zeta_\pm^L}{1-\sum_{i=(+,-)} \zeta_i^L}\right) \pm \phi^L$$

$$\zeta_\pm^L = c_\pm^L m_\pm \qquad c_\pm^L = c_\pm(x = L)$$

$$\phi^L = \phi(x = L)$$

## Concentration distribution in a size-modified case (channel)

## **Usage**
### **Running the main script**
```sh
julia --project=. ./src/smpnp.jl input.yml
```
**Arguments:**
- `--paramaters` or `-p`: Path to the input data file.
- `--help` or `-h`: Display help information

### **Example**
```sh
python src/main.py --input data/example.csv --output results/processed.csv
```

---

## **Results & Output**
This code generates:
- **Processed data:** Stored in `results/`
- **Figures & plots:** Saved as `.png` files
- **Log files:** Stored in `logs/` (if applicable)

---

## **Citation**
If you use this code in your research, please cite:

```bibtex
@article{YourName2024,
  author  = {Your Name and Co-Authors},
  title   = {Your Paper Title},
  journal = {Journal Name},
  year    = {2024},
  volume  = {XX},
  pages   = {XX-XX},
  doi     = {DOI HERE}
}
```
or:
> Your Name et al., *Your Paper Title*, in **Journal Name**, 2024. DOI: [Insert DOI].

---

## **License**
This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

## **Contact**
For questions or collaboration inquiries, please contact:

* Niklaus M. Leuenberger,  
Stanford University, [Department of Energy Science & Engineering](https://ese.stanford.edu/)  
Green Earth Sciences Bldg. Rm 151  
367 Panama Street, Stanford, CA 94305, U.S.A  
[niklausl@stanford.edu](mailto:niklausl@stanford.edu)
 
## **References**

[^bazant_lecture_notes]: Lecture Notes [MIT-10.626 Lecture notes](https://ocw.mit.edu/courses/10-626-electrochemical-energy-systems-spring-2014/34aaca3a97887695dd295db7cc0fa3c0_MIT10_626S14_S11lec24.pdf)


