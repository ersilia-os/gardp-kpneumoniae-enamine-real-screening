import os
import h5py
import time
import httplib2
import pandas as pd
import numpy as np
from tqdm import tqdm
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build

GDRIVE_FOLDER_ID = "1FBELagBf9hlKVgvkaZ8YF60jKRAmsHPo"
ARRAY_DTYPE = "int8"


def _get_data_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_X.npz")


def _get_identifiers_file(dir_path, chunk_name):
    return os.path.join(dir_path, f"{chunk_name}_SMILES_IDs.tsv.zip")


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
    service_file = os.path.abspath(os.path.join(root, "..", "..", "config", "service.json"))
    folder_id = GDRIVE_FOLDER_ID
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
    t0 = time.time()
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
    t1 = time.time()
    print(f"Download completed in {t1 - t0:.2f} seconds.")


def convert_to_h5(dir_path, chunk_name, batch_size=100_000):
    t0 = time.time()
    print(f"Converting chunk {chunk_name} to H5 format...")
    data_file = _get_data_file(dir_path, chunk_name)
    identifiers_file = _get_identifiers_file(dir_path, chunk_name)
    h5_file = _get_h5_file(dir_path, chunk_name)
    identifiers = pd.read_csv(identifiers_file, compression="zip", delimiter="\t")
    print(identifiers.head())
    X = np.load(data_file)["X"]
    n_rows, n_cols = X.shape
    print(f"Identifiers rows: {len(identifiers)}, Data shape: {n_rows}, {n_cols}")
    assert len(identifiers) == n_rows, "Row mismatch between X and identifiers"
    with h5py.File(h5_file, "w") as f:
        dset_X = f.create_dataset("values", shape=(n_rows, n_cols), dtype=ARRAY_DTYPE, chunks=(batch_size, n_cols), compression="gzip")
        f.create_dataset("key", data=identifiers["id"].tolist(), dtype=h5py.string_dtype())
        f.create_dataset("input", data=identifiers["smiles"].tolist(), dtype=h5py.string_dtype())
        for start in tqdm(range(0, n_rows, batch_size), desc="Converting to H5"):
            end = min(start + batch_size, n_rows)
            X_batch = X[start:end].astype(ARRAY_DTYPE)
            dset_X[start:end] = X_batch
            del X_batch
    t1 = time.time()
    print(f"H5 file saved to {h5_file}")
    print(f"Conversion completed in {t1 - t0:.2f} seconds.")
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