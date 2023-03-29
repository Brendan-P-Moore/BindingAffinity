
import pandas as pd
import numpy as np
import Bio

import rdkit
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import IPythonConsole
from rdkit import DataStructs
from rdkit.DataStructs import ConvertToNumpyArray

import requests as r
from Bio import SeqIO
from io import StringIO
from chembl_webresource_client.new_client import new_client

import multiprocessing
from multiprocessing import Pool, cpu_count
import itertools

import sklearn
from sklearn.model_selection import train_test_split, GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr

import gc
from functools import reduce

import peptides
import propy
from propy import PyPro

import joblib
from joblib import dump, load
import os

from autogluon.features.generators import AutoMLPipelineFeatureGenerator
from autogluon.tabular import FeatureMetadata
from autogluon.tabular import TabularDataset, TabularPredictor




def uniprot_to_fasta(df):
    uniprot_ids = df['target_id'].drop_duplicates().tolist()
    baseUrl="https://rest.uniprot.org/uniprotkb/"
    seq = []
    for id in uniprot_ids:
        currentUrl = baseUrl + id + ".fasta"
        response = r.get(currentUrl)
        Seq = StringIO(response.text)
        pSeq = list(SeqIO.parse(Seq,'fasta'))
        if not pSeq:  # Check if pSeq is empty
            seq.append("")  # Append an empty string if no sequence was found
            print(f"Warning: No sequence found for Uniprot ID {id}.")
        else:
            seq_records = []
            for seq_record in pSeq:
                seq_records.append(str(seq_record.seq))
                break
            seq.append(seq_records[0])
    df['sequence'] = df['target_id'].map(dict(zip(uniprot_ids,seq)))
    df.dropna(inplace=True)
    df.reset_index(drop=True,inplace=True)
    return df 


def process_batch(batch):
    molecules = new_client.molecule.filter(molecule_chembl_id__in=batch).only('molecule_chembl_id','molecule_structures')
    smiles=[]
    chembl_ids_with_smiles = []
    for molecule in molecules:
        if molecule['molecule_structures'] is not None and 'canonical_smiles' in molecule['molecule_structures']:
            smiles.append(molecule['molecule_structures']['canonical_smiles'])
            chembl_ids_with_smiles.append(molecule['molecule_chembl_id'])
    return dict(zip(chembl_ids_with_smiles, smiles))


def chembl_to_smiles(df):
    chembl_ids = df['compound_id'].unique().tolist()
    batch_size = 500
    chembl_id_batches = [chembl_ids[i:i+batch_size] for i in range(0, len(chembl_ids), batch_size)]
    num_processes = cpu_count()
    with Pool(num_processes) as pool:
        results = pool.map(process_batch, chembl_id_batches)
    merged_results = reduce(lambda x, y: {**x, **y}, results)
    df['compound_smiles'] = df['compound_id'].map(merged_results)
    df.dropna(inplace=True)
    df.reset_index(drop=True,inplace=True)
    return df


def smi_to_mol(smi):
            return Chem.MolFromSmiles(smi)
    
def smiles_to_rdkit_mol(df):            
    with Pool(cpu_count()) as p:
        mols = p.map(smi_to_mol, df['compound_smiles'])
    mols = [mol for mol in mols if mol is not None]
    return mols


def fingerprint_generation(df,ncounts): # returns a df with molecular fingerprints, ncounts is the number of count vectors.
    mols = smiles_to_rdkit_mol(df)
    # fingerprint generation
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=ncounts) # morgan300 fingerprint
    ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=ncounts) # TopologicalTorsion300 Fingerprint
    mfp_matrix = np.zeros((len(mols), ncounts), dtype=int)
    tt_matrix = np.zeros((len(mols), ncounts), dtype=int)
    mfp_fingerprints = [mfpgen.GetCountFingerprint(mol) for mol in mols]
    tt_fingerprints = [ttgen.GetCountFingerprint(mol) for mol in mols]
    for i, fingerprint in enumerate(mfp_fingerprints):
        ConvertToNumpyArray(fingerprint, mfp_matrix[i])
    mfp_df = pd.DataFrame(mfp_matrix, columns=[f"morg_count_fp_{i}" for i in range(ncounts)]).astype('int8')
    for i, fingerprint in enumerate(tt_fingerprints):
        ConvertToNumpyArray(fingerprint, tt_matrix[i])
    tt_df = pd.DataFrame(tt_matrix, columns=[f"tt_fp_{i}" for i in range(ncounts)]).astype('int8')
    merged_df = pd.concat([df,mfp_df,tt_df],axis=1)
    return merged_df


def calculate_descriptors(unique_seqs, calculation_function):
    unique_descriptors = {}
    for seq in unique_seqs:
        unique_descriptors[seq] = calculation_function(seq)
    descriptor_names = set().union(*unique_descriptors.values())
    descriptor_dict = {}
    for name in descriptor_names:
        descriptor_dict[name] = {}
        for seq, desc in unique_descriptors.items():
            descriptor_dict[name][seq] = desc.get(name, None)
    return descriptor_dict

def protein_feature_generation(df):
    unique_seqs = set(df['sequence'])
    
    # Calculate descriptors using peptides library with full code, different syntax required.
    unique_descriptors = {}
    for seq in unique_seqs:
        unique_descriptors[seq] = peptides.Peptide(seq).descriptors()

    descriptor_names = set().union(*unique_descriptors.values())
    descriptor_dict = {}
    for name in descriptor_names:
        descriptor_dict[name] = {}
        for seq, desc in unique_descriptors.items():
            descriptor_dict[name][seq] = desc.get(name, None)
    peptide_df = pd.DataFrame.from_dict(descriptor_dict)
    peptide_df['sequence'] = peptide_df.index
    
    # Calculate descriptors using propy.CTD library
    unique_descriptors = calculate_descriptors(unique_seqs, propy.CTD.CalculateCTD)
    ctd_df = pd.DataFrame.from_dict(unique_descriptors)
    ctd_df['sequence'] = ctd_df.index
    
    # Calculate descriptors using propy.Autocorrelation library
    unique_descriptors = calculate_descriptors(unique_seqs, propy.Autocorrelation.CalculateGearyAutoTotal)
    autocorr_df = pd.DataFrame.from_dict(unique_descriptors)
    autocorr_df['sequence'] = autocorr_df.index
    
    # Calculate descriptors using propy.AAComposition library
    unique_descriptors = calculate_descriptors(unique_seqs, propy.AAComposition.CalculateAAComposition)
    aa_df = pd.DataFrame.from_dict(unique_descriptors)
    aa_df['sequence'] = aa_df.index
    
    # Merge all dataframes into one
    df = pd.merge(df, peptide_df, on='sequence')
    df = pd.merge(df, ctd_df, on='sequence')
    df = pd.merge(df, autocorr_df, on='sequence')
    df = pd.merge(df, aa_df, on='sequence')
    
    df['sequence_length'] = df['sequence'].str.len()
    
    return df.copy()


def autogluon_predict(df):
    exclude_autogluon= ['compound_id','target_id','group','sequence','compound_smiles']
    features_autogluon = [c for c in df.columns if c not in exclude_autogluon]
    autogluon_folds=10
    model = TabularPredictor.load(f"autogluon_models/fold0/")
    pred_autogluon = model.predict(df[features_autogluon])
    for f in range(1, autogluon_folds):
        model = TabularPredictor.load(
            f"autogluon_models/fold{f}/"
        )
        pred_autogluon += model.predict(df[features_autogluon])
    pred_autogluon /= autogluon_folds
    df['predicted_pKd'] = pred_autogluon
    return df
    
# For a dataframe containing a column "compound_id" in chembl format, and a column "target_id" in UniprotID format
def predict_from_chembl_uniprot(df):
    df = uniprot_to_fasta(df)
    df = chembl_to_smiles(df)
    df = fingerprint_generation(df,300)
    df = protein_feature_generation(df)
    df = autogluon_predict(df)
    return df

def predict_from_fasta_smiles(df):
    df = fingerprint_generation(df,300)
    df = protein_feature_generation(df)
    df = autogluon_predict(df)
    return df