import os
import numpy as np
import joblib
import json
import time
import h5py
import pandas as pd
from tqdm import tqdm


root = os.path.dirname(os.path.abspath(__file__))


STRAINS = ["abaumannii", "kpneumoniae_nctc_13438", "kpneumoniae_atcc_43816"]
CUTOFFS = ["bin", "perc01", "perc05", "perc1"]
MODEL_TYPES = ["full", "attractive"]

GARDP_ENDPOINTS = []
for strain in STRAINS:
    for cutoff in CUTOFFS:
        for model_type in MODEL_TYPES:
            model_name = f"{strain}_{cutoff}_{model_type}"
            GARDP_ENDPOINTS += [model_name]

OTHER_ACTIVITY_ENDPOINTS = [
    "stokes_abaumannii",
    "stokes_ecoli",
    "mole_gn",
    "gneprop_tolc",
]

TRANSPORT_ENDPOINTS = [
    "entry_rules",
    "gn_permeability_proxy",
]

TOXICITY_ENDPOINTS = [
    "cytotoxicity"
]

ATTRACTIVENESS_ENDPOINTS = [
    "antibioticdb",
    "collins_abx",
    "abx_resemblance",
    "frequent_hitters",
    "pains"
]


KEEP_ENDPOINTS = GARDP_ENDPOINTS + OTHER_ACTIVITY_ENDPOINTS + TRANSPORT_ENDPOINTS

DISCARD_ENDPOINTS = TOXICITY_ENDPOINTS + ATTRACTIVENESS_ENDPOINTS


class LightDecisionModel(object):

    def __init__(self, name, model, cutoff, keep):
        print(f"Initializing LightDecisionModel: {name}, cutoff={cutoff}, keep={keep}")
        self.name = name
        self.model = model
        self.cutoff = cutoff
        self.keep = keep
        self.chunk_size = 1000000

    def screen(self, h5_file, idxs):
        print(f"Screening with model: {self.name} on H5 file {h5_file}")
        if type(idxs) is list:
            idxs = np.array(idxs, dtype="int")
        chunk_size = self.chunk_size
        y_hat = np.zeros(len(idxs), dtype=np.float16)
        with h5py.File(h5_file, "r") as f:
            X = f["values"]
            for i in tqdm(range(0, len(idxs), chunk_size)):
                chunk_idxs = idxs[i:i+chunk_size]
                X_chunk = X[chunk_idxs,:]
                y_hat[i:i+len(chunk_idxs)] = self.model.predict_proba(X_chunk)[:, 1]
        print(f"Median predicted score: {np.median(y_hat)}, cutoff: {self.cutoff}")
        if self.keep:
            y_keep = y_hat >= self.cutoff
        else:
            y_keep = y_hat <= self.cutoff
        print(f"Number of compounds kept: {np.sum(y_keep)} out of {len(y_keep)}")
        return y_keep


class LightScreener(object):

    def __init__(self):
        self.models_dir = os.path.abspath(os.path.join(root, "..", "..", "data", "endpoints"))
        self.models = {}
        for endpoint in GARDP_ENDPOINTS + OTHER_ACTIVITY_ENDPOINTS + TRANSPORT_ENDPOINTS + TOXICITY_ENDPOINTS + ATTRACTIVENESS_ENDPOINTS:
            model = self._load_model(endpoint)
            self.models[endpoint] = model
        
    def _load_model(self, name):
        model_json = os.path.join(self.models_dir, f"{name}.json")
        model_joblib = os.path.join(self.models_dir, f"{name}.joblib")
        with open(model_json, "r") as f:
            data = json.load(f)
            cutoff = data["cutoff"]
        model = joblib.load(model_joblib)
        if name in KEEP_ENDPOINTS:
            keep = True
        elif name in DISCARD_ENDPOINTS:
            keep = False
        else:
            raise Exception(f"Unknown endpoint: {name}")
        return LightDecisionModel(name, model, cutoff, keep)
    
    def _get_identifiers(self, h5_file):
        with h5py.File(h5_file, "r") as f:
            key_list = [x.decode("utf-8") for x in f["key"][:]]
            smiles_list = [x.decode("utf-8") for x in f["input"][:]]
        data = {
            "key": key_list,
            "smiles": smiles_list
        }
        return data

    def _screen(self, h5_file):
        print("###########################################################")
        print("Starting GARDP LightScreener...")
        t0 = time.time()

        idxs = np.arange(len(X), dtype="int")
        print("- Initial number of compounds:", len(X))

        print("Cytotoxicity screening...")
        model = self.models["cytotoxicity"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after cytotoxicity screening:", len(idxs))

        print("PAINS filter...")
        model = self.models["pains"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after PAINS filter:", len(idxs))

        print("Frequent hitters...")
        model = self.models["frequent_hitters"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after frequent hitters filter:", len(idxs))

        print("Collins similarity...")
        model = self.models["collins_abx"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after Collins similarity filter:", len(idxs))

        print("AntibioticDB similarity...")
        model = self.models["antibioticdb"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after AntibioticDB similarity filter:", len(idxs))

        print("Antibiotic resemblance...")
        model = self.models["abx_resemblance"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after antibiotic resemblance filter:", len(idxs))

        print("Abaumannii GARDP screening with bin cutoff...")
        model = self.models["abaumannii_bin_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_bin_full = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii bin full filter:", len(idxs_abaumannii_bin_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with perc01 cutoff...")
        model = self.models["abaumannii_perc01_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_perc01_full = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii perc01 full filter:", len(idxs_abaumannii_perc01_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with bin cutoff (attractive)...")
        model = self.models["abaumannii_bin_attractive"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_bin_attractive = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii bin attractive filter:", len(idxs_abaumannii_bin_attractive))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with perc01 cutoff (attractive)...")
        model = self.models["abaumannii_perc01_attractive"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_perc01_attractive = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii perc01 attractive filter:", len(idxs_abaumannii_perc01_attractive))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with bin cutoff...")
        model = self.models["kpneumoniae_nctc_13438_bin_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_bin_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 bin full filter:", len(idxs_kpneumoniae_nctc_13438_bin_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc01 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc01_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_perc01_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 perc01 full filter:", len(idxs_kpneumoniae_nctc_13438_perc01_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))
        
        print("Kpneumoniae NCTC13438 GARDP screening with bin cutoff (attractive)...")
        model = self.models["kpneumoniae_nctc_13438_bin_attractive"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_bin_attractive = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 bin attractive filter:", len(idxs_kpneumoniae_nctc_13438_bin_attractive))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc01 cutoff (attractive)...")
        model = self.models["kpneumoniae_nctc_13438_perc01_attractive"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_perc01_attractive = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 perc01 attractive filter:", len(idxs_kpneumoniae_nctc_13438_perc01_attractive))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("#######################################################")

        print("Abaumannii GARDP screening with perc05 cutoff...")
        model = self.models["abaumannii_perc05_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_perc05_full = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii perc05 full filter:", len(idxs_abaumannii_perc05_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...:", len(idxs))

        print("Abaumannii GARDP screening with perc1 cutoff...")
        model = self.models["abaumannii_perc1_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_abaumannii_perc1_full = idxs[keep_mask]
        print("- Number of compounds kept for the A.baumannii perc1 full filter:", len(idxs_abaumannii_perc1_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...:", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc05 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc05_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_perc05_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 perc05 full filter:", len(idxs_kpneumoniae_nctc_13438_perc05_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...:", len(idxs))
        
        print("Kpneumoniae NCTC13438 GARDP screening with perc1 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc1_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_nctc_13438_perc1_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae NCTC13438 perc1 full filter:", len(idxs_kpneumoniae_nctc_13438_perc1_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...:", len(idxs))

        print("Kpneumoniae ATCC43816 GARDP screening with bin cutoff...")
        model = self.models["kpneumoniae_atcc_43816_bin_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_atcc_43816_bin_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae ATCC43816 bin full filter:", len(idxs_kpneumoniae_atcc_43816_bin_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae ATCC43816 GARDP screening with perc01 cutoff...")
        model = self.models["kpneumoniae_atcc_43816_perc01_full"]
        keep_mask = model.screen(h5_file, idxs)
        idxs_kpneumoniae_atcc_43816_perc01_full = idxs[keep_mask]
        print("- Number of compounds kept for the K.pneumoniae ATCC43816 perc01 full filter:", len(idxs_kpneumoniae_atcc_43816_perc01_full))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        idxs = []
        idxs += list(idxs_abaumannii_perc05_full)
        idxs += list(idxs_abaumannii_perc1_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc05_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc1_full)
        idxs += list(idxs_kpneumoniae_atcc_43816_bin_full)
        idxs += list(idxs_kpneumoniae_atcc_43816_perc01_full)
        idxs = sorted(idxs)
        idxs = np.array(idxs, dtype="int")
        print("Compounds left for testing after GARDP screening:", len(idxs))

        print("#######################################################")

        print("Other activity filters...")

        print("Stokes E.coli filter...")
        model = self.models["stokes_ecoli"]
        keep_mask = model.screen(X)
        idxs_stokes_ecoli = idxs[keep_mask]
        print("- Number of compounds kept for the Stokes E.coli filter:", len(idxs_stokes_ecoli))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Stokes A.baumannii filter...")
        model = self.models["stokes_abaumannii"]
        keep_mask = model.screen(X)
        idxs_stokes_abaumannii = idxs[keep_mask]
        print("- Number of compounds kept for the Stokes A.baumannii filter:", len(idxs_stokes_abaumannii))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("Mol-E Gram-negative filter...")
        model = self.models["mole_gn"]
        keep_mask = model.screen(X)
        idxs_mole_gn = idxs[keep_mask]
        print("- Number of compounds kept for the Mole GN filter:", len(idxs_mole_gn))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        print("GNEProp TolC filter...")
        model = self.models["gneprop_tolc"]
        keep_mask = model.screen(X)
        idxs_gneprop_tolc = idxs[keep_mask]
        print("- Number of compounds kept for the GNEProp TolC filter:", len(idxs_gneprop_tolc))
        idxs = idxs[~keep_mask]
        print("- Number of compounds left for testing...", len(idxs))

        idxs = []
        idxs += list(idxs_stokes_ecoli)
        idxs += list(idxs_stokes_abaumannii)
        idxs += list(idxs_mole_gn)
        idxs += list(idxs_gneprop_tolc)
        idxs = sorted(idxs)
        idxs = np.array(idxs)

        print("####################################################") 
        model = self.models["entry_rules"]
        keep_mask = model.screen(h5_file, idxs)
        idxs = idxs[keep_mask]
        print("- Number of compounds after entry rules filter:", len(idxs))

        model = self.models["gn_permeability_proxy"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        print("- Number of compounds after GN permeability proxy filter:", len(idxs))

        idxs = list(idxs)
        idxs += list(idxs_abaumannii_bin_full)
        idxs += list(idxs_abaumannii_perc01_full)
        idxs += list(idxs_abaumannii_bin_attractive)
        idxs += list(idxs_abaumannii_perc01_attractive)
        idxs += list(idxs_kpneumoniae_nctc_13438_bin_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc01_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_bin_attractive)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc01_attractive)
        idxs = sorted(idxs)
        idxs = np.array(idxs, dtype="int")

        print("####################################################")
        print("Final number of compounds selected:", len(idxs))
        t1 = time.time()
        print(f"Total screening time: {t1 - t0:.2f} seconds")
        
        return idxs
    
    def screen(self, h5_input, csv_output):
        identifiers_data = self._get_identifiers(h5_input)
        idxs = self._screen(h5_input)
        smiles_list = [identifiers_data["smiles"][i] for i in idxs]
        key_list = [identifiers_data["key"][i] for i in idxs]
        df = pd.DataFrame({
            "identifier": key_list,
            "smiles": smiles_list
        })
        df.to_csv(csv_output, index=False)
        print(f"Results saved to {csv_output}")

