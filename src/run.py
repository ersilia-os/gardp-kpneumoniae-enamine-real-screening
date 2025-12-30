import argparse

from core.core import pipeline
from core.utils import download_data, copy_data, convert_to_h5, clean_data


def main():

    parser = argparse.ArgumentParser(description="Run the GARDP LightScreener on a specified data chunk.")
    parser.add_argument("--chunk_name", type=str, required=True, help="Name of the data chunk to process.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output CSV file.")
    parser.add_argument("--from_dir", type=str, default=None, required=False, help="Directory to fetch fingerprints from. If provided, it skips downloading.")
    parser.add_argument("--skip_clean", action='store_true', help="If set, skips the cleaning step after processing.")
    args = parser.parse_args()

    chunk_name = args.chunk_name
    output_dir = args.output_dir
    from_dir = args.from_dir
    skip_clean = args.skip_clean

    if from_dir is None:
        download_data(output_dir, chunk_name)
    else:
        copy_data(from_dir, output_dir, chunk_name)
    convert_to_h5(output_dir, chunk_name)
    pipeline(output_dir, chunk_name)
    
    if not skip_clean:
        clean_data(output_dir, chunk_name)


if __name__ == "__main__":
    main()