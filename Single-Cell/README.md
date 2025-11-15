## R Environment Notes

- Packages from **CRAN** are installed using `install.packages()`.
- Packages from **Bioconductor** are installed using `BiocManager::install()`.
- Packages hosted on **GitHub** are installed using `remotes::install_github()`.

No additional environment configuration or customization is required beyond these standard installation methods.  
**RStudio** and **RStudio Server** are used as the primary IDEs.

## Python Environment Notes

- The environment is primarily managed using **miniconda**, which provides the core dependencies.
- Most packages are installed through **pip**, along with any remaining dependencies not covered by conda.
- An `environment.yaml` file can be provided upon request to reproduce the exact same environment specification.
- **Jupyter Notebook** and **Spyder** are commonly used as IDEs.

## Code Numbering Convention

The numerical prefixes in file names indicate the intended execution order.  
Within a pipeline, scripts with **smaller numbers should be executed before those with larger numbers** to ensure proper workflow sequencing.
