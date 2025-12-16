import os
import tempfile
import argparse

from core.core import LightScreener
from core.utils import download_data, convert_to_h5, clean_data


def main():

    parser = argparse.ArgumentParser(description="Run the GARDP LightScreener on a specified data chunk.")
    parser.add_argument("--chunk_name", type=str, help="Name of the data chunk to process.")
    parser.add_argument("--output_dir", type=str, help="Directory to save the output CSV file.")
    parser.add_argument("--gdrive_api_key", type=str, help="Google Drive API key for downloading data.")
    args = parser.parse_args()

    chunk_name = args.chunk_name
    output_dir = args.output_dir
    gdrive_api_key = args.gdrive_api_key

    tmp_dir = tempfile.mkdtemp(prefix="ersilia-")
    download_data(chunk_name, tmp_dir, gdrive_api_key)
    h5_file = convert_to_h5(chunk_name, tmp_dir)
    screener = LightScreener()
    screener.screen(h5_input=h5_file, csv_output=os.path.join(output_dir, f"{chunk_name}_hits.csv"))
    clean_data(chunk_name, tmp_dir)


if __name__ == "__main__":
    main()