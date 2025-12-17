import os
import h5py
import httplib2
import pandas as pd
import numpy as np
from scipy.sparse import load_npz
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build


def _get_data_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_X.npz")


def _get_identifiers_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_SMILES_IDs.csv.zip")


def _get_h5_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}.h5")


def check_exists(dir_path, chunk_name):
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    if os.path.exists(data_file) and os.path.exists(identifiers_file):
        return True
    else:
        return False
    
def download_file(outfile):
    root = os.path.dirname(os.path.abspath(__file__))
    file = os.path.basename(outfile)
    service_file = os.path.abspath(os.path.join(root, "..", "..", "data", "service", "service.json"))
    folder_id = "19x9yAUySBXgrHBE3gjGjHsomcLQuLQmL"
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.readonly"])
    service = build("drive", "v3", credentials=creds)
    httplib2.DEFAULT_TIMEOUT = 600

    query = f"name='{file}' and '{folder_id}' in parents and trashed=false"
    results = service.files().list(q=query, fields="files(id)", supportsAllDrives=True, includeItemsFromAllDrives=True).execute()
    files = results.get("files", [])
    if not files:
        raise FileNotFoundError(f"'{file}' not found! Consider checking available chunks in data/chunks/chunks.csv")
    if len(files) > 1:
        raise RuntimeError(
            f"Multiple files named '{file}' are found...")

    file_id = files[0]["id"]
    request = service.files().get_media(fileId=file_id, supportsAllDrives=True)
    with open(outfile, "wb") as f:
        f.write(request.execute())


def download_data(dir_path, chunk_name):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    if os.path.exists(data_file):
        os.remove(data_file)
    if os.path.exists(identifiers_file):
        os.remove(identifiers_file)
    print(f"Downloading identifiers_file to {identifiers_file}...")
    download_file(identifiers_file)
    print(f"Downloading data_file to {data_file}...")
    download_file(data_file)


def convert_to_h5(dir_path, chunk_name, batch_size=100_000):
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    h5_file = _get_h5_file(dir_path, chunk_name)
    X_sparse = load_npz(data_file)   # csr_matrix
    n_rows, n_cols = X_sparse.shape
    identifiers = pd.read_csv(identifiers_file, compression="zip")
    print(len(identifiers), n_rows)
    assert len(identifiers) == n_rows, "Row mismatch between X and identifiers"
    with h5py.File(h5_file, "w") as f:
        dset_X = f.create_dataset("values", shape=(n_rows, n_cols), dtype="int8", chunks=(batch_size, n_cols), compression="gzip")
        f.create_dataset("key", data=np.array(identifiers["IDS"]).astype("S"))
        f.create_dataset("input", data=np.array(identifiers["SMILES"]).astype("S"))
        for start in range(0, n_rows, batch_size):
            end = min(start + batch_size, n_rows)
            print(f"Writing rows {start}:{end}")
            X_batch = X_sparse[start:end].toarray().astype("int8")
            dset_X[start:end] = X_batch
            del X_batch
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