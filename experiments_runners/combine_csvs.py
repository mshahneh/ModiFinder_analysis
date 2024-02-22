import pandas as pd
import os
import sys
import argparse


def main(direcotry_path, output_path):
    # get all csv files in the directory
    csv_files = [
        os.path.join(direcotry_path, f)
        for f in os.listdir(direcotry_path)
        if f.endswith(".csv")
    ]

    # read all csv files
    dfs = [pd.read_csv(f) for f in csv_files]

    # combine all csv files into one
    df = pd.concat(dfs)

    # write to csv file
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="combine all csv files in a directory into one csv file"
    )
    parser.add_argument(
        "--directory_path", type=str, help="path to the directory containing csv files"
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="scores.csv",
        help="path to the output file",
    )

    args = parser.parse_args()

    main(args.directory_path, args.output_path)
