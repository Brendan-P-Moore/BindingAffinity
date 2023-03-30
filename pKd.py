import pandas as pd
import numpy as np
import requests as r
import gc
import multiprocessing
from multiprocessing import Pool, cpu_count
import joblib
from joblib import dump, load
import os
import tqdm
from tqdm import tqdm
import itertools

import rdkit
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import IPythonConsole
from rdkit import DataStructs
from rdkit.DataStructs import ConvertToNumpyArray

import Bio
from Bio import SeqIO
from io import StringIO
from chembl_webresource_client.new_client import new_client
import peptides
import propy
from propy import PyPro

import sklearn
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr

from autogluon.features.generators import AutoMLPipelineFeatureGenerator
from autogluon.tabular import FeatureMetadata
from autogluon.tabular import TabularDataset, TabularPredictor


def drug_target_pairs(sequences_path, smiles_path):
    fasta_df = pd.read_csv(os.path.abspath(sequences_path), names=["sequence"])

    smiles_df = pd.read_csv(os.path.abspath(smiles_path), names=["compound_smiles"])

    combinations = list(
        itertools.product(fasta_df["sequence"], smiles_df["compound_smiles"])
    )

    df = pd.DataFrame(combinations, columns=["sequence", "compound_smiles"])
    return df


def uniprot_to_fasta(df):
    uniprot_ids = df["target_id"].drop_duplicates().tolist()
    baseUrl = "https://rest.uniprot.org/uniprotkb/"
    seq = []
    for id in tqdm(uniprot_ids, desc="Converting uniprot IDs to protein sequences"):
        currentUrl = baseUrl + id + ".fasta"
        response = r.get(currentUrl)
        Seq = StringIO(response.text)
        pSeq = list(SeqIO.parse(Seq, "fasta"))
        if not pSeq:  # Check if pSeq is empty
            seq.append("")  # Append an empty string if no sequence was found
            print(f"Warning: No sequence found for Uniprot ID {id}.")
        else:
            seq_records = []
            for seq_record in pSeq:
                seq_records.append(str(seq_record.seq))
                break
            seq.append(seq_records[0])
    df["sequence"] = df["target_id"].map(dict(zip(uniprot_ids, seq)))
    df.dropna(inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def process_batch(batch):
    molecules = new_client.molecule.filter(molecule_chembl_id__in=batch).only(
        "molecule_chembl_id", "molecule_structures"
    )
    smiles = []
    chembl_ids_with_smiles = []
    for molecule in molecules:
        if (
            molecule["molecule_structures"] is not None
            and "canonical_smiles" in molecule["molecule_structures"]
        ):
            smiles.append(molecule["molecule_structures"]["canonical_smiles"])
            chembl_ids_with_smiles.append(molecule["molecule_chembl_id"])
    return dict(zip(chembl_ids_with_smiles, smiles))


def chembl_to_smiles(df):
    chembl_ids = df["compound_id"].unique().tolist()
    batch_size = 500
    chembl_id_batches = [
        chembl_ids[i : i + batch_size] for i in range(0, len(chembl_ids), batch_size)
    ]
    num_processes = cpu_count()
    with Pool(num_processes) as pool:
        results = []
        for result in tqdm(
            pool.imap_unordered(process_batch, chembl_id_batches),
            total=len(chembl_id_batches),
        ):
            results.append(result)
    merged_results = {}
    for d in results:
        merged_results.update(d)
    df["compound_smiles"] = df["compound_id"].map(merged_results)
    df.dropna(inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def smi_to_mol(smi):
    return Chem.MolFromSmiles(smi)


def smiles_to_rdkit_mol(df):
    with Pool(cpu_count()) as p:
        mols = p.map(smi_to_mol, df["compound_smiles"])
    mols = [mol for mol in mols if mol is not None]
    return mols


def fingerprint_generation(
    df, ncounts
):  # returns a df with molecular fingerprints, ncounts is the number of count vectors.
    mols = smiles_to_rdkit_mol(df)
    # fingerprint generation
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=ncounts
    )  # morgan300 fingerprint
    ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(
        fpSize=ncounts
    )  # TopologicalTorsion300 Fingerprint
    mfp_matrix = np.zeros((len(mols), ncounts), dtype=int)
    tt_matrix = np.zeros((len(mols), ncounts), dtype=int)
    mfp_fingerprints = [
        mfpgen.GetCountFingerprint(mol)
        for mol in tqdm(mols, desc="Generating Morgan fingerprints", unit=" molecule")
    ]
    tt_fingerprints = [
        ttgen.GetCountFingerprint(mol)
        for mol in tqdm(
            mols, desc="Generating Topological Torsion fingerprints", unit=" molecule"
        )
    ]

    for i, fingerprint in enumerate(mfp_fingerprints):
        ConvertToNumpyArray(fingerprint, mfp_matrix[i])
    mfp_df = pd.DataFrame(
        mfp_matrix, columns=[f"morg_count_fp_{i}" for i in range(ncounts)]
    ).astype("int8")

    for i, fingerprint in enumerate(tt_fingerprints):
        ConvertToNumpyArray(fingerprint, tt_matrix[i])
    tt_df = pd.DataFrame(
        tt_matrix, columns=[f"tt_fp_{i}" for i in range(ncounts)]
    ).astype("int8")

    merged_df = pd.concat([df, mfp_df, tt_df], axis=1)
    return merged_df


def calculate_descriptors(unique_seqs, calculation_function):
    unique_descriptors = {}
    for seq in tqdm(unique_seqs, desc="Calculating propy package descriptors"):
        unique_descriptors[seq] = calculation_function(seq)
    descriptor_names = set().union(*unique_descriptors.values())
    descriptor_dict = {}
    for name in descriptor_names:
        descriptor_dict[name] = {}
        for seq, desc in unique_descriptors.items():
            descriptor_dict[name][seq] = desc.get(name, None)
    return descriptor_dict


def protein_feature_generation(df):
    unique_seqs = set(df["sequence"])

    # Calculate descriptors using peptides library with full code, different syntax required.
    unique_descriptors = {}
    for seq in tqdm(unique_seqs, desc="Calculating all peptide package descriptors"):
        unique_descriptors[seq] = peptides.Peptide(seq).descriptors()

    descriptor_names = set().union(*unique_descriptors.values())
    descriptor_dict = {}
    for name in descriptor_names:
        descriptor_dict[name] = {}
        for seq, desc in unique_descriptors.items():
            descriptor_dict[name][seq] = desc.get(name, None)
    peptide_df = pd.DataFrame.from_dict(descriptor_dict)
    peptide_df["sequence"] = peptide_df.index

    # Calculate descriptors using propy.CTD library
    unique_descriptors = calculate_descriptors(unique_seqs, propy.CTD.CalculateCTD)

    ctd_df = pd.DataFrame.from_dict(unique_descriptors)
    ctd_df["sequence"] = ctd_df.index

    # Calculate descriptors using propy.Autocorrelation library
    unique_descriptors = calculate_descriptors(
        unique_seqs, propy.Autocorrelation.CalculateGearyAutoTotal
    )
    autocorr_df = pd.DataFrame.from_dict(unique_descriptors)
    autocorr_df["sequence"] = autocorr_df.index

    # Calculate descriptors using propy.AAComposition library
    unique_descriptors = calculate_descriptors(
        unique_seqs, propy.AAComposition.CalculateAAComposition
    )
    aa_df = pd.DataFrame.from_dict(unique_descriptors)
    aa_df["sequence"] = aa_df.index

    # Merge all dataframes into one
    df = pd.merge(df, peptide_df, on="sequence")
    df = pd.merge(df, ctd_df, on="sequence")
    df = pd.merge(df, autocorr_df, on="sequence")
    df = pd.merge(df, aa_df, on="sequence")

    df["sequence_length"] = df["sequence"].str.len()

    return df


def autogluon_predict(df):
    exclude_autogluon = [
        "compound_id",
        "target_id",
        "group",
        "sequence",
        "compound_smiles",
    ]
    features_autogluon = [c for c in df.columns if c not in exclude_autogluon]
    autogluon_models = 10
    model = TabularPredictor.load(f"autogluon_models/fold0/")
    pred_autogluon = model.predict(df[features_autogluon])
    with tqdm(total=autogluon_models - 1, desc="Autogluon Model Predictions") as pbar:
        for f in range(1, autogluon_models):
            model = TabularPredictor.load(f"autogluon_models/fold{f}/")
            pred_autogluon += model.predict(df[features_autogluon])
            pbar.update(1)
    pred_autogluon /= autogluon_models
    df["predicted_pKd"] = pred_autogluon
    return df


# For a dataframe containing a column "compound_id" in chembl format, and a column "target_id" in UniprotID format
def pKd_from_chembl_uniprot(df):
    df = uniprot_to_fasta(df)
    df = chembl_to_smiles(df)
    df = fingerprint_generation(df, 300)
    df = protein_feature_generation(df)
    df = autogluon_predict(df)
    return df


# For a dataframe containing a column "compound_smiles" in smiles format, and a column "sequence" of protein sequences for all drug-target pairs
def pKd_from_sequence_smiles(df):
    df = fingerprint_generation(df, 300)
    df = protein_feature_generation(df)
    df = autogluon_predict(df)
    return df


# For a text file of compounds in smiles notation, and a separate text file of protein sequences. Each should be a single column with no header.
def pKd_from_txt(sequences_path, smiles_path):
    df = drug_target_pairs(sequences_path, smiles_path)
    df = fingerprint_generation(df, 300)
    df = protein_feature_generation(df)
    df = autogluon_predict(df)
    return df
