# QuantumWalk-LongCovid

create_SIP_network

A Python pipeline to build a SIP (Systematic Interaction Profile) network using DIP, DEP, and STRING resources.

Defaults: If any of the above keys are missing, safe defaults will be applied (e.g. half of your CPU cores for num_processes, 10000 for chunk_size, etc.).

🛠️ Dependencies
	•	Python 3.7+
	•	pandas
	•	networkx

Install via:

pip install pandas networkx

pandas>=1.0
networkx>=2.5

🚀 Usage
	1.	Clone or download the repository.
	2.	Place your config.json in the project root (or update paths inside it).


python create_SIP_network.py

All outputs (graphs in Pickle format, logs, etc.) will appear under the directory specified by result_dir.

📝 License

Distributed under the MIT License. See LICENSE for details.
