import argparse

from core.core import pipeline
from core.utils import download_data, convert_to_h5, clean_data


def main():

    parser = argparse.ArgumentParser(description="Run the GARDP LightScreener on a specified data chunk.")
    parser.add_argument("--chunk_name", type=str, required=True, help="Name of the data chunk to process.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output CSV file.")
    args = parser.parse_args()

    chunk_name = args.chunk_name
    output_dir = args.output_dir

    download_data(output_dir, chunk_name)
    convert_to_h5(output_dir, chunk_name)
    pipeline(output_dir, chunk_name)
    clean_data(output_dir, chunk_name)


if __name__ == "__main__":
    main()