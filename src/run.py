import os
import argparse

from core.core import LightScreener
from core.utils import download_data, convert_to_h5, clean_data


def main():

    parser = argparse.ArgumentParser(description="Run the GARDP LightScreener on a specified data chunk.")
    parser.add_argument("--chunk_name", type=str, required=True, help="Name of the data chunk to process.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output CSV file.")
    args = parser.parse_args()

    chunk_name = args.chunk_name
    output_dir = args.output_dir

    #download_data(output_dir, chunk_name)
    #h5_file = convert_to_h5(output_dir, chunk_name)
    h5_file = os.path.join(output_dir, f"{chunk_name}.h5") # TODO remove when debugged
    screener = LightScreener()
    screener.screen(h5_input=h5_file, csv_output=os.path.join(output_dir, f"{chunk_name}_hits.csv"))
    #clean_data(output_dir, chunk_name)


if __name__ == "__main__":
    main()