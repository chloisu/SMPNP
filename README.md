# **SMPNP - a size modified Poisson Nernst Planck simulation tool**

<!-- [![DOI](https://zenodo.org/badge/DOI.svg)](https://doi.org/your-doi)  -->
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)  
[![Python Version](https://img.shields.io/badge/python-3.10.16%2B-blue)](#prerequisites)  

## **\ud83d\udcda Overview**
This repository contains the source code for the research project **"[Your Research Paper Title]"**, which was published in **[Journal/Conference Name]**. The code is designed to **[briefly describe functionality, e.g., numerical simulations, data processing, machine learning, etc.]**.

\ud83d\udcdd **Paper DOI**: [Insert DOI]  
\ud83d\udcc1 **Data DOI (if available)**: [Insert DOI for dataset]  
\ud83d\udd17 **Project Website (if applicable)**: [Insert link]  

---

## **\ud83d\udcc2 Repository Structure**
```
├── src/                 # Main source code
│   ├── main.py          # Entry point script
│   ├── utils.py         # Utility functions
│   ├── analysis.py      # Data analysis
├── notebooks/           # Jupyter notebooks for analysis
├── data/                # Sample datasets (not included in repo)
├── results/             # Output results
├── requirements.txt     # Python dependencies
├── environment.yml      # Conda environment file (if applicable)
├── pyenv_setup.sh       # Pyenv setup script
├── Dockerfile           # Docker container setup (if applicable)
├── README.md            # This file
└── LICENSE              # License file
```

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

## **\u25b6\ufe0f Usage**
### **Running the main script**
```sh
python src/main.py --input data/sample.csv --output results/output.csv
```
**Arguments:**
- `--input`: Path to the input data file.
- `--output`: Path to save the processed results.

### **Example**
```sh
python src/main.py --input data/example.csv --output results/processed.csv
```

---

## **\ud83d\udcca Results & Output**
This code generates:
- **Processed data:** Stored in `results/`
- **Figures & plots:** Saved as `.png` files
- **Log files:** Stored in `logs/` (if applicable)

---

## **\ud83d\udcdd Citation**
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


