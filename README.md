# QuantumWalk-LongCovid


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This project is a bioinformatics pipeline that utilizes a **Quantum Walk** algorithm to identify genes associated with Long COVID. It integrates Protein-Protein Interaction (PPI) data to construct a SIP network, then runs a quantum walk simulation on it to predict and prioritize genes highly associated with Long COVID.

The complete pipeline consists of the following main stages:

1.  **Data Integration and Network Construction**: Integrates PPI data from multiple sources such as DIP, DEP, and STRING to generate a SIP network.
2.  **Quantum Walk Simulation**: Executes a quantum walk algorithm on the constructed network to calculate the importance of each node (gene).
3.  **Results Analysis and Gene Identification**: Analyzes the results of the quantum walk to produce a list of genes predicted to be most relevant to Long COVID.

### 🛠️ Dependencies

To run this pipeline, you will need the following libraries:
* Python 3.7+
* `pandas`
* `networkx`
* `numpy`
* `scipy`
* `hiperwalk`

You can install the necessary libraries using `pip`:
```bash
pip install pandas networkx numpy scipy hiperwalk
```

### 🚀 Usage

**1. Clone Repository**

First, clone this repository to your local machine.
```bash
git clone https://github.com/Namshik-Han-Lab/QuantumWalk-LongCovid.git
cd QuantumWalk-LongCovid
```

**2. Prepare Data**

Place all required input data files in the `data/` directory. The paths and filenames specified in the `config.json` file must match the actual files.

Please download the PPI link data from STRING DB, which excluded from the Git repository due to large size.

https://stringdb-downloads.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz

**3. Configuration**

Open the `config.json` file in the project's root directory to configure the pipeline's execution environment. An example `config.json` is provided below:

**Key Configuration Details:**
* `dip_file`, `dep_file`, `string_info_file`, `string_links_file`: Paths to the input data files used for network construction.
* `result_dir`: The directory where all outputs (network files, logs, gene lists, etc.) will be saved.
* `num_processes`: The number of CPU cores to use for data processing.
* `combined_score`: The minimum score (threshold) to use when filtering interactions from the STRING database.
* `chunk_size`: The number of rows to read into memory at a time when processing large files.
* (If there are parameters related to the quantum walk, additional explanations should be added here.)

**4. Run Pipeline**

The entire pipeline is executed by running several scripts in sequence. (The commands below should be adjusted according to the actual script names in the project.)

**Step 1: Create SIP Network**
```bash
python create_SIP_network_with_covid19_edges.py
```
This script uses the data specified in `config.json` to create a `sip_network.pkl` file in the `result/` directory.

**Step 2: Run Quantum Walk**
```bash
python network_walking.py
```
This script loads the generated SIP network, performs the quantum walk simulation, and saves `COVID_RWR_result.pkl` and `COVID_coined_result.pkl` in the `result/` directory.

**Step 3: Generate mixed Result**
```bash
python generate_mixed_result.py
```
This script loads the generated network walking results, calculates the quantum walk mixing time, and saves `quantum_walk_result.csv` in the `result/` directory.

### 📂 Output

Once the pipeline execution is complete, the following files will be generated in your `result_dir`:

* **`sip_network.pkl`**: The integrated protein-protein interaction network file. It can be loaded using the `networkx` library.
* **`COVID_RWR_result.pkl`**: Random walk results as Dataframe. 
* **`COVID_coined_result.pkl`**: Quantum walk results as Dataframe.
* **`quantum_walk_result.csv`**: A CSV file containing the quantum walk simulation results for each gene (e.g., final probability, rank).
* **`pipeline.log`**: A log file that records the entire execution process of the pipeline, including warnings and errors.

### 📝 License

This project is distributed under the MIT License. For more details, see the `LICENSE` file.

### 📜 Citation

If you find this work useful, we would appreciate it if you cite the relevant publication. (Add citation information here once the paper is published.)
