# Multiple Sequence Alignment (MSA)

## 🧬 Project Overview
This project explores **extensions of the Needleman-Wunsch algorithm** to support multiple sequence alignment (MSA). The main goal is to study and compare different ways of adapting the classical dynamic programming approach to align more than two biological sequences.

The project implements:
- **Progressive alignment** of multiple sequences
- **Merging of pre-aligned sequence groups**
- **Exact multidimensional alignment**
- Use of **Carillo-Lipman bounds** to optimize search space

Despite algorithmic improvements, **exact multidimensional alignment remains computationally expensive** (∼ N^K), making it unsuitable for large datasets. Carillo-Lipman bounds **can reduce the number of computed cells**, although the effect depends on input sequences.

---

## ⚙️ Alignment Principles

All algorithms implemented in this project share the same biological alignment logic:

- **Gap opening cost > gap extension cost** to favor fewer, longer gaps (biologically plausible).
- Support for both **nucleotide (DNA/RNA)** and **protein** sequences.
  
### Scoring Parameters

| Mode              | Match | Mismatch | Gap Opening | Gap Extension | Matrix     |
|-------------------|-------|----------|-------------|---------------|------------|
| DNA/RNA           | +1    | -1       | -10         | -2            | Simple     |
| Protein (BLOSUM)  | --    | --       | -10         | -2            | BLOSUM62   |

To use BLOSUM62, set `blosum_m=True` in relevant functions.

---

## 🧠 Implemented Algorithms

### ✔️ Pairwise Needleman-Wunsch
- Classic global alignment between two sequences.

### ✔️ Progressive N-sequence alignment
- Sequences are aligned sequentially:
  - First, \( S_0 \) and \( S_1 \)
  - Then \( S_2 \) to the previous alignment, and so on

### ✔️ N+M Group Alignment
- Two pre-aligned sequence sets are merged.
- Enables modular construction of large alignments.

### ✔️ Exact Multidimensional Alignment
- Constructs a full \( N^K \) dynamic programming matrix.
- Aligns all sequences simultaneously (not progressively).
- Limited to small \( K \) due to exponential complexity.

---

## 📊 Evaluation and Comparison

- **Correctness validation**: pairwise alignments are compared with **Biopython**'s reference output.
- **Quality metric**: **Sum-of-pairs (SP) score**, commonly used for MSA evaluation.

---

## 🔍 Matrix and Alignment Visualization

Visual tools are provided to explore algorithm behavior and assist debugging:

- **2D and 3D matrix plots**
- **Traceback path visualizations**
- **Alignment reconstruction display**

For \( K = 3 \), **heatmap slices** of the 3D DP matrix can be generated to observe scoring patterns and path trends.

---

## ⚡ Carillo-Lipman Bounds

The **Carillo-Lipman bounding method** is used to reduce computational load in the multidimensional case. It defines **valid index ranges** for each sequence based on pairwise alignments and skips computation outside these bounds.

- **Precomputes feasible regions**
- **Skips invalid transitions**
- **Assigns \( -\infty \)** to unreachable cells

This speeds up the DP matrix filling phase, especially when alignment paths are well-constrained.

---

## ⚠️ Limitations

Despite various optimizations, the implemented algorithms have several inherent limitations:

- **Exponential complexity**: The multidimensional Needleman-Wunsch algorithm has time and space complexity of \( O(N^K) \), making it impractical for aligning more than 4–5 sequences.
- **Scalability issues**: Progressive and group alignment strategies scale better, but they rely on alignment order and may lead to suboptimal global results.
- **CL bounds effectiveness varies**: Carillo-Lipman constraints reduce computation only when sequences share high similarity or structure. In other cases, bounds may have little impact.
- **Limited visualization**: Matrix visualization is currently only feasible for up to 3 sequences due to dimensionality constraints.
- **No parallelization**: Current implementation runs on a single CPU thread and is not optimized for performance on large datasets.

These limitations highlight the challenges of exact MSA methods and motivate the use of heuristics in practical applications.

---

## 📂 Project Structure
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

---

## 🚀 Installation
Clone the repository and install dependencies:
```bash
git clone https://github.com/akochergina/multiple-seq-alignment.git
cd multiple-seq-alignment
```

---

## ⚙️ Usage
To explore the implemented methods and their applications, run the Jupyter Notebook located in the `notebooks/` directory:
```bash
jupyter notebook notebooks/using_example.ipynb
```
This notebook provides a **step-by-step demonstration** of the implemented functions.