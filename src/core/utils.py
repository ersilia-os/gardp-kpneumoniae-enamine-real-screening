import os
import h5py
import pandas as pd
import numpy as np


def _get_data_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_X.npz")


def _get_identifiers_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_SMILES_IDs.csv.zst")


def _get_h5_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}.h5")


def check_exists(dir_path, chunk_name):
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    if os.path.exists(data_file) and os.path.exists(identifiers_file):
        return True
    else:
        return False


def download_data(dir_path, chunk_name, gdrive_api_key):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    if os.path.exists(data_file):
        os.remove(data_file)
    if os.path.exists(identifiers_file):
        os.remove(identifiers_file)
    # TODO Arnau: implement actual download logic here


def convert_to_h5(dir_path, chunk_name):
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    h5_file = _get_h5_file(dir_path, chunk_name)
    X = np.load(data_file)['X']
    identifiers = pd.read_csv(identifiers_file)
    with h5py.File(h5_file, 'w') as f:
        f.create_dataset('values', data=X.astype("int8"))
        f.create_dataset("key", data=np.array(identifiers["IDS"]).astype("S"))
        f.create_dataset("input", data=np.array(identifiers["SMILES"]).astype("S"))
    return h5_file


def clean_data(dir_path, chunk_name):
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    h5_file = _get_h5_file(dir_path, chunk_name)
    if os.path.exists(data_file):
        os.remove(data_file)
    if os.path.exists(identifiers_file):
        os.remove(identifiers_file)
    if os.path.exists(h5_file):
        os.remove(h5_file)