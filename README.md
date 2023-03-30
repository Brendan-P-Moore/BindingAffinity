# BindingAffinity
BindingAffinity utilizes proteochemometric (PCM) modeling to predict the binding affinity for drug-protein pairs. In this context, the binding affinity is the -log(Kd) (pKd), where Kd is the dissociation constant for a drug-protein pair. The model was trained using molecular fingerprints combined with physicochemical protein descriptors, and was tuned and tested on chemically dissimilar datasets (1). The task of predicting chemically dissimilar drug-protein pairs is much more challenging and realistic, as it requires the model to generalize to novel chemical space where new drug-protein interactions can be discovered (1). This predictor improves on the baseline put forward by Guvenilir and Doğan. Details can be found in the description section.


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
git clone https://github.com/Brendan-P-Moore/BindingAffinity

```
2. Set the current working directory to the BindingAffinity folder.

```
cd BindingAffinity

```
3. Install the required packages listed in the dependencies section and requirements file.

```
pip install -r requirements.txt

```
4. Import the pKd python file, and the predict function.

```
import pKd
from pKd import predict_pKd_from_txt

```
# Predict




```
predict_pKd_from_txt(sequences_path,smiles_path)

```


# Description




# References

1) Atas Guvenilir, H., Doğan, T. How to approach machine learning-based prediction of drug/compound–target interactions. J Cheminform 15, 16 (2023). https://doi.org/10.1186/s13321-023-00689-w
