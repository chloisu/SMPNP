# **SMPNP - a size modified Poisson Nernst Planck simulation tool**

<!-- [![DOI](https://zenodo.org/badge/DOI.svg)](https://doi.org/your-doi)  -->
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Python Version](https://img.shields.io/badge/python-3.10.16%2B-blue)](#prerequisites)  

## **Overview**
This repository contains the source code for the research project **"[Your Research Paper Title]"**, which was published in **[Journal/Conference Name]**. The code is designed to **[briefly describe functionality, e.g., numerical simulations, data processing, machine learning, etc.]**.

\ud83d\udcdd **Paper DOI**: [Insert DOI]  
\ud83d\udcc1 **Data DOI (if available)**: [Insert DOI for dataset]  
\ud83d\udd17 **Project Website (if applicable)**: [Insert link]  

---

## **Installation**
### **Prerequisites**
Before running the code, ensure you have the following installed:
- **Python \u2265 3.10** (recommended: `pyenv` or `conda`)
- **Git** for cloning the repository
- **Required libraries** (see `requirements.txt`)

### **2\ufe0f\u20e3 Setup Instructions**
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
$$\psi_e \quad \text{electrical potential}$$

### **Non-dimensionalization**
The variables in governing equations are dimensionless. They are made non-dimensional by the following scalings wich can also be found in[^bazant_lecture_notes]

$\tilde{x} = \frac{x}{L}, \quad \tilde{c_i} = \frac{c_i}{c_{i,ref}}, \quad \tilde{\varepsilon} = \frac{\varepsilon}{\varepsilon_{ref}},\quad \tilde{D_i} = \frac{D_i}{D_{i,ref}}, \quad \tilde{t} =\frac{t}{\left(L^2/D\right)},\quad \tilde{\phi} = \frac{e\phi}{k_BT}, \tilde{\nabla} = L\nabla$

In the following, we drop the tilda notation for non-dimensional variables and assume that all variables are non-dimensional.
### **Governing equations**
The dimensionless governing PNP equations then read as
$$
    \frac{\partial c^-}{\partial t} = \nabla \cdot \left(M_- c_-\nabla \mu_-(\mathbf{c})\right)\\
    \frac{\partial \mathbf{c}_+}{\partial t} = \nabla \cdot \left(M_+ \mathbf{c}_+ \nabla \mu_+(\mathbf{c})\right)\\
    \lambda_D \nabla^2 \phi = -\frac{1}{\varepsilon}(c_+-c_-)
$$

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
This project is licensed under the **[License Name]**. See the [LICENSE](LICENSE) file for details.

## **Contact**
For questions or collaboration inquiries, please contact:

* Niklaus M. Leuenberger,  
Stanford University, [Department of Energy Science & Engineering](https://ese.stanford.edu/)  
Green Earth Sciences Bldg. Rm 151  
367 Panama Street, Stanford, CA 94305, U.S.A  
[niklausl@stanford.edu](mailto:niklausl@stanford.edu)
 
## **References**

[^bazant_lecture_notes]: Lecture Notes [MIT-10.626 Lecture notes](https://ocw.mit.edu/courses/10-626-electrochemical-energy-systems-spring-2014/34aaca3a97887695dd295db7cc0fa3c0_MIT10_626S14_S11lec24.pdf)


