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

BindingAffinity is a machine learning model for the prediction of drug-protein binding affinities. It was trained on the newly published benchmark datasets produced by Guvenilir and Doğan (1). Guvenilir and Doğan generated 3 different splits of train/test data for 10 different protein families (1):

1.	Randomly split data, which produces highly overoptimistic testing results.
2.	A dissimilar compound split, in which only dissimilar drug molecules were included in the testing data.
3.	A fully dissimilar split, where both proteins and compounds in the test set were not chemically similar to the training data.

Splits were prepared in this way for each of the 10 protein families. The fully dissimilar split proved much more challenging for predictors, as it requires more generalizable predictions. This testing scenario is also more realistic for cases where novel drug-protein pairs are desired.
This package combines all 10 protein family training sets for one larger training set spanning a greater range of chemical space. It then trains 10 different autogluon models, each one tuned on a different withheld family of proteins. The trained autogluon models are ensembles of PyTorch tabular networks and lightgbm models. Because different protein families have dissimilar drug-protein interactions, this withheld tuning data is suitable as a cross-validation dataset for the chemically dissimilar testing data. The figure below summarizes this methodology.

![alt text](https://github.com/Brendan-P-Moore/BindingAffinity/blob/main/training_scheme_BindingAffinity.png?raw=true)

This predictor shows an improvement over the baseline random forest models trained by Doğan and Guvenilir, with an improvement of 0.062 RMSE, and 0.17 higher spearman correlation coefficient compared with the best single model of Guvenilir and Doğan. This shows that with chemically dissimilar tuning data, machine learning models can better predict novel drug-protein binding affinities, despite the use of easy to generate molecular fingerprints and physicochemical descriptors. The table below shows a detailed comparison of the best overall baseline model from Guvenilir and Doğan, a baseline random forest model trained on the molecular and protein descriptors used in this work, and the final predictor. The enzyme family tuned ensemble model greatly outperforms the baseline models in terms of spearman rank correlation on all enzyme families. The root mean squared error was also lower overall for the tuned model, although epigenetic regulators and ion-channels error was lower in the baseline model.

| Spearman   Correlation Coefficient | Overall | Proteases | Epigenetic regulators | Hydrolases | Ion-channel | membrane receptors | other | oxidoreductases | transcription-factors | transferases | transporters |
|------------------------------------|---------|-----------|-----------------------|------------|-------------|--------------------|-------|-----------------|-----------------------|--------------|--------------|
| Guvenilir   Random Forest          | 0.302   | 0.413     | 0.445                 | 0.280      | 0.153       | 0.341              | 0.349 | 0.079           | 0.266                 | 0.439        | 0.249        |
| Random Forest Baseline             | 0.393   | 0.352     | 0.492                 | 0.451      | 0.231       | 0.299              | 0.329 | 0.145           | 0.339                 | 0.352        | 0.225        |
| Autogluon Tuned                    | 0.476   | 0.449     | 0.465                 | 0.472      | 0.291       | 0.425              | 0.437 | 0.256           | 0.320                 | 0.533        | 0.372        |
|                                    |         |           |                       |            |             |                    |       |                 |                       |              |              |
| Root Mean Squared Error            | Overall | Proteases | Epigenetic regulators | Hydrolases | Ion-channel | membrane receptors | other | oxidoreductases | transcription-factors | transferases | transporters |
| Guvenilir   Random Forest          | 1.230   | 1.207     | 0.975                 | 1.676      | 1.144       | 1.229              | 1.347 | 1.271           | 1.150                 | 1.132        | 1.173        |
| Random Forest Baseline             | 1.251   | 1.297     | 0.953                 | 1.368      | 1.110       | 1.313              | 1.252 | 1.208           | 1.150                 | 1.247        | 1.196        |
| Autogluon Tuned                    | 1.168   | 1.171     | 1.041                 | 1.304      | 1.155       | 1.206              | 1.208 | 1.115           | 1.098                 | 1.119        | 1.100        |

Another key issue with proteochemometric (PCM) based binding affinity predictors is the tendency to ignore protein features and over-rely on molecular features. Since lightgbm models played an important role in the final predictor ensemble, the feature importance of lightgbm models trained using the scheme shown above were investigated. The figure below shows the top 20 feature importance of lightgbm models tuned on (A): ion-channels, (B): membrane receptors, (C): transcription factors, and (D): transferases. We can see that while a few molecular fingerprint features were typically most important, the top 20 features include a good mix of protein features and molecular features. Features labelled morg_count are from the morgan300 count vector, and tt_count are from the topological torsion 300 count vector. Features starting with GearyAuto are geary autocorrelation features generated suing propy3, and the remaining features are either peptides/propy3 descriptors, or in the case of single letters amino acid distribution features.

![alt text](https://github.com/Brendan-P-Moore/BindingAffinity/blob/main/feature_importance_lightgbm_BindingAffinity.png?raw=true)

# References

1) Atas Guvenilir, H., Doğan, T. How to approach machine learning-based prediction of drug/compound–target interactions. J Cheminform 15, 16 (2023). https://doi.org/10.1186/s13321-023-00689-w

2) https://github.com/hevalatas/ProtBENCH
