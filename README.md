# Drug-Protein-pKd-Predictor
Predicts the dissociation coefficients between drug-protein (target) pairs using molecular fingerprints and physicochemical protein properties. The molecular fingerprints used are the morgan 300 countvector fingerprint, and 300 topological torsion countvector fingerprint.


# Dependencies
* python >= 3.9
* autogluon==0.7.0
* chembl_webresource_client == 0.10.8
* biopython==1.81
* propy3 == 1.1.1
* peptides == 0.3.1
* rdkit==2022.9.5

# Installation Guide

1. Clone the GitHub repository.

```
git clone https://github.com/Brendan-P-Moore/Drug-Protein-pKd-Predictor

```
2. Set the current working directory to the Drug-Protein-pKd-Predictor folder.

```
cd Drug-Protein-pKd-Predictor

```
3. Install the required packages listed in the dependencies section and requirements file.

```
pip install -r requirements.txt

```
4. Import the AMPTransformer python file, and the predict function.

```
import pKd
from pKd import predict_pKd_from_txt

```
# Predict




```
predict_pKd_from_txt(sequences_path,smiles_path)

```
