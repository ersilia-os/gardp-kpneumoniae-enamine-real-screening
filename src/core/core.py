import os
import numpy as np
import joblib
import json
import time
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

    def screen(self, X):
        print(f"Screening with model: {self.name}")
        chunk_size = self.chunk_size
        num_chunks = int(np.ceil(X.shape[0] / chunk_size))
        y_hat = np.zeros(X.shape[0])
        for i in tqdm(range(num_chunks)):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, X.shape[0])
            y_hat[start_idx:end_idx] = self.model.predict_proba(X[start_idx:end_idx])[:, 1]
        y_hat = self.model.predict_proba(X)[:,1]
        print(f"Median predicted score: {np.median(y_hat)}, cutoff: {self.cutoff}")
        if self.keep:
            y_keep = y_hat >= self.cutoff
        else:
            y_keep = y_hat <= self.cutoff
        return y_keep
    
    def screen_h5(self, h5_file):



class LightScreener(object):

    def __init__(self):
        self.models_dir = os.path.join(root, "..", "..", "data", "endpoints")
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
        identifiers = None
        return identifiers

    def _screen(self, h5_file):

        X_original = X.copy()

        t0 = time.time()

        idxs = np.arange(len(X))
        print("Initial number of compounds:", len(X))

        print("Cytotoxicity screening...")
        model = self.models["cytotoxicity"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after cytotoxicity screening:", len(idxs))

        print("PAINS filter...")
        model = self.models["pains"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after PAINS filter:", len(idxs))

        print("Frequent hitters...")
        model = self.models["frequent_hitters"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after frequent hitters filter:", len(idxs))

        print("Collins similarity...")
        model = self.models["collins_abx"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after Collins similarity filter:", len(idxs))

        print("AntibioticDB similarity...")
        model = self.models["antibioticdb"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after AntibioticDB similarity filter:", len(idxs))

        print("Antibiotic resemblance...")
        model = self.models["abx_resemblance"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after antibiotic resemblance filter:", len(idxs))

        print("Abaumannii GARDP screening with bin cutoff...")
        model = self.models["abaumannii_bin_full"]
        keep_mask = model.screen(X)
        idxs_abaumannii_bin_full = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii bin full filter:", len(idxs_abaumannii_bin_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with perc01 cutoff...")
        model = self.models["abaumannii_perc01_full"]
        keep_mask = model.screen(X)
        idxs_abaumannii_perc01_full = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii perc01 full filter:", len(idxs_abaumannii_perc01_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with bin cutoff (attractive)...")
        model = self.models["abaumannii_bin_attractive"]
        keep_mask = model.screen(X)
        idxs_abaumannii_bin_attractive = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii bin attractive filter:", len(idxs_abaumannii_bin_attractive))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Abaumannii GARDP screening with perc01 cutoff (attractive)...")
        model = self.models["abaumannii_perc01_attractive"]
        keep_mask = model.screen(X)
        idxs_abaumannii_perc01_attractive = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii perc01 attractive filter:", len(idxs_abaumannii_perc01_attractive))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with bin cutoff...")
        model = self.models["kpneumoniae_nctc_13438_bin_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_bin_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 bin full filter:", len(idxs_kpneumoniae_nctc_13438_bin_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc01 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc01_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_perc01_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 perc01 full filter:", len(idxs_kpneumoniae_nctc_13438_perc01_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))
        
        print("Kpneumoniae NCTC13438 GARDP screening with bin cutoff (attractive)...")
        model = self.models["kpneumoniae_nctc_13438_bin_attractive"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_bin_attractive = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 bin attractive filter:", len(idxs_kpneumoniae_nctc_13438_bin_attractive))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc01 cutoff (attractive)...")
        model = self.models["kpneumoniae_nctc_13438_perc01_attractive"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_perc01_attractive = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 perc01 attractive filter:", len(idxs_kpneumoniae_nctc_13438_perc01_attractive))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("#######################################################")

        print("Abaumannii GARDP screening with perc05 cutoff...")
        model = self.models["abaumannii_perc05_full"]
        keep_mask = model.screen(X)
        idxs_abaumannii_perc05_full = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii perc05 full filter:", len(idxs_abaumannii_perc05_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...:", len(idxs))

        print("Abaumannii GARDP screening with perc1 cutoff...")
        model = self.models["abaumannii_perc1_full"]
        keep_mask = model.screen(X)
        idxs_abaumannii_perc1_full = idxs[keep_mask]
        print("Number of compounds kept for the A.baumannii perc1 full filter:", len(idxs_abaumannii_perc1_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...:", len(idxs))

        print("Kpneumoniae NCTC13438 GARDP screening with perc05 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc05_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_perc05_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 perc05 full filter:", len(idxs_kpneumoniae_nctc_13438_perc05_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...:", len(idxs))
        
        print("Kpneumoniae NCTC13438 GARDP screening with perc1 cutoff...")
        model = self.models["kpneumoniae_nctc_13438_perc1_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_nctc_13438_perc1_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae NCTC13438 perc1 full filter:", len(idxs_kpneumoniae_nctc_13438_perc1_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...:", len(idxs))

        print("Kpneumoniae ATCC43816 GARDP screening with bin cutoff...")
        model = self.models["kpneumoniae_atcc_43816_bin_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_atcc_43816_bin_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae ATCC43816 bin full filter:", len(idxs_kpneumoniae_atcc_43816_bin_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        print("Kpneumoniae ATCC43816 GARDP screening with perc01 cutoff...")
        model = self.models["kpneumoniae_atcc_43816_perc01_full"]
        keep_mask = model.screen(X)
        idxs_kpneumoniae_atcc_43816_perc01_full = idxs[keep_mask]
        print("Number of compounds kept for the K.pneumoniae ATCC43816 perc01 full filter:", len(idxs_kpneumoniae_atcc_43816_perc01_full))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        idxs = []
        idxs += list(idxs_abaumannii_perc05_full)
        idxs += list(idxs_abaumannii_perc1_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc05_full)
        idxs += list(idxs_kpneumoniae_nctc_13438_perc1_full)
        idxs += list(idxs_kpneumoniae_atcc_43816_bin_full)
        idxs += list(idxs_kpneumoniae_atcc_43816_perc01_full)
        idxs = sorted(idxs)
        idxs = np.array(idxs)
        X = X_original[idxs,:]
        print("Compounds left for testing...")

        print("#######################################################")

        model = self.models["stokes_ecoli"]
        keep_mask = model.screen(X)
        idxs_stokes_ecoli = idxs[keep_mask]
        print("Number of compounds kept for the Stokes E.coli filter:", len(idxs_stokes_ecoli))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        model = self.models["stokes_abaumannii"]
        keep_mask = model.screen(X)
        idxs_stokes_abaumannii = idxs[keep_mask]
        print("Number of compounds kept for the Stokes A.baumannii filter:", len(idxs_stokes_abaumannii))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        model = self.models["mole_gn"]
        keep_mask = model.screen(X)
        idxs_mole_gn = idxs[keep_mask]
        print("Number of compounds kept for the Mole GN filter:", len(idxs_mole_gn))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        model = self.models["gneprop_tolc"]
        keep_mask = model.screen(X)
        idxs_gneprop_tolc = idxs[keep_mask]
        print("Number of compounds kept for the GNEProp TolC filter:", len(idxs_gneprop_tolc))
        idxs = idxs[~keep_mask]
        X = X[~keep_mask]
        print("Number of compounds left for testing...", len(idxs))

        idxs = []
        idxs += list(idxs_stokes_ecoli)
        idxs += list(idxs_stokes_abaumannii)
        idxs += list(idxs_mole_gn)
        idxs += list(idxs_gneprop_tolc)
        idxs = sorted(idxs)
        idxs = np.array(idxs)
        X = X_original[idxs,:]

        print("####################################################") 
        model = self.models["entry_rules"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after entry rules filter:", len(idxs))

        model = self.models["gn_permeability_proxy"]
        keep_mask = model.screen(X)
        idxs = idxs[keep_mask]
        X = X[keep_mask]
        print("Number of compounds after GN permeability proxy filter:", len(idxs))

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

        print("Final number of compounds selected:", len(idxs))
        t1 = time.time()
        print(f"Total screening time: {t1 - t0:.2f} seconds")
        
        return idxs
    
    def screen(self, h5_input, csv_output):
        identifiers = self._get_identifiers(h5_input)
        idxs = self._screen(h5_input)
        df = pd.DataFrame(identifiers[idxs], columns=[])
        df.to_csv(csv_output, index=False)

