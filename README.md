# Multiple Sequence Alignment (MSA)

## Project Overview
This project focuses on implementing and optimizing the **Needleman-Wunsch algorithm** for multiple sequence alignment (MSA). Various enhancements and extensions of the classical method are explored, including **Carillo-Lipman bounds** for reducing computational complexity.

Despite optimizations, the algorithm **does not scale well** due to its high computational complexity. Additionally, **Carillo-Lipman bounds do not always significantly reduce the search space**, depending on the dataset. The project primarily focuses on **developing and evaluating different improvements and extensions** of the Needleman-Wunsch algorithm.

## Project Structure
```
/project_root/
    ├── balibase/                       # Real data from BALIBASE for testing 
    ├── modules/                        # Python modules with isolated functions
    │   ├── needleman_wunsch/           # Needleman-Wunsch functions
    │   │   ├── nw_core.py              # Core Needleman-Wunsch implementations
    │   │   ├── nw_cost.py              # Cost functions for Needleman-Wunsch algorithm
    │   │   ├── nw_visualization.py     # Visualization of matrices and alignment results
    │   ├── carillo_lipman_bounds.py    # Carrillo-Lipman bounds for optimization
    │   ├── result_evaluation.py        # Functions to compare different alignments 
    ├── notebooks/ 
    │   ├── using_example.ipynb         # Example Jupyter Notebook demonstrating key functions
    ├── README.md                       # Project documentation
```

## Installation
Clone the repository and install dependencies:
```bash
git clone https://github.com/akochergina/multiple-seq-alignment.git
cd multiple-seq-alignment
```

## Usage
To explore the implemented methods and their applications, run the Jupyter Notebook located in the `notebooks/` directory:
```bash
jupyter notebook notebooks/using_example.ipynb
```
This notebook provides a **step-by-step demonstration** of the implemented functions.

## Features
- **Multiple sequence alignment** using an extended **Needleman-Wunsch algorithm**.
- **Carillo-Lipman bounds optimization**, which can reduce the computational space in some cases.
- **Backtracking and traceback visualization** for reconstructing optimal alignments.
- **Evaluation functions** to compare different alignment results.

## Limitations
- The **algorithm has exponential complexity** in the number of sequences (∼ N^K), making it **unsuitable for large-scale alignments**.
- **Carrillo-Lipman bounds do not always significantly reduce computational space**, as their effectiveness depends on the sequences being aligned.
