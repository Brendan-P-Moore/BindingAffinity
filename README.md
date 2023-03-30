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

# Prediction Guide

1. Import the pKd python file. To predict from two text files, one a file with a single protein sequence on each line, and one a file with a single smiles string on each line, import the "predict_pKd_from_txt" function. Example sequences and smiles text files can be found in the "Examples" folder in this repository.

```
import pKd
from pKd import pKd_from_txt

```

2. The prediction can be made for every combination of peptide and compound found in the two text files by specifying the absolute path to the files. To avoid memory issues, limit the file sizes to less than 50 proteins, and less than 1000 compounds per prediction. This prediction returns a dataframe with all features used to predict the binding affinity, and a column titled "predicted_pKd" with the predicted -log(Kd) value.

```
prediction_dataframe = pKd_from_txt(sequences_path,smiles_path)

predicted_pKd = prediction_dataframe['predicted_pKd']

```

Using the example files provided in this repository:

```
prediction_dataframe = pKd_from_txt("examples/sequence.txt","examples/smiles.txt")

predicted_pKd = prediction_dataframe['predicted_pKd']

```

3. Alternatively, the prediction can be made directly from a dataframe of drug-protein pairs. It is highly recommended that the proteins are represented by their sequences, and the drugs by their smiles strings for faster predictions. The sequences column must be named "sequence", and the smiles string column must be named "compound_smiles". df in the following example is the sequences/smiles dataframe.

```
import pKd
from pKd import pKd_from_sequence_smiles

prediction_dataframe = pKd_from_sequence_smiles(df)

predicted_pKd = prediction_dataframe['predicted_pKd']

```

Binding affinities can also be predicted from a dataframe of drugs in chembl id format and proteins in uniprot id format. Below is a worked example using the test data which is found in this format. This prediction will be much slower, as the chembl ids need to be converted to SMILES, and the uniprot ids need to be converted to sequences.

```
import pKd
from pKd import pKd_from_chembl_uniprot

df = pd.read_csv('test_data/test_uniprot_chembl.csv')

prediction_dataframe = pKd_from_chembl_uniprot(df)

predicted_pKd = prediction_dataframe['predicted_pKd']

```

# Description



# References

1) Atas Guvenilir, H., Doğan, T. How to approach machine learning-based prediction of drug/compound–target interactions. J Cheminform 15, 16 (2023). https://doi.org/10.1186/s13321-023-00689-w

2) https://github.com/hevalatas/ProtBENCH
