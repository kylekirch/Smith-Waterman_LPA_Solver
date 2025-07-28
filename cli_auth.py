"""
Filename: cli_auth.py
Project: Smith-Waterman Local Pairwise Alignment (LPA)
Description: Handles command-line interface authentication and validates
              input/output/scoring matrix file paths.
Author: Kyle Kirchgessner
Date: 2025-07-28
Version: 1.0
Dependencies: argparse, os, sys
"""


import argparse
import os
import sys


def cli_auth():
    '''
    Description:
    input file handling

    IN:
    command line argument

    OUT:
    path to an input file (-i)
    path to an output file (-o)
    path to a .mtx file (-s)
    '''

    def validate_input_file(path):
        # verify that the input file exists
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                f"Input file does not exist: {path}")
        return path

    def validate_score_file(path):
        # verify that the input file exists
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                f"Matrix file does not exist: {path}")
        return path

    def validate_output_path(path):
        # verify that the output file exists
        dir_path = os.path.dirname(
            path) or '.'  # falls back to the current directory as a root
        if not os.path.isdir(dir_path):
            raise argparse.ArgumentTypeError(
                f"Directory for output file does not exist: {dir_path}")
        return path

    parser = argparse.ArgumentParser(
        description="Validate input and output file paths.")

    parser.add_argument(
        "-i",
        type=validate_input_file,
        required=True,
        help="Path to input file (must exist)."
    )

    parser.add_argument(
        "-o",
        type=validate_output_path,
        required=True,
        help="Path to output file (directory must exist)."
    )

    parser.add_argument(
        "-s",
        type=validate_score_file,
        required=True,
        help="Path to Sore Matrix file (must exist)."
    )
    return (parser.parse_args())


def main():
    args = cli_auth()

    input_file_path = args.i
    output_file_path = args.o
    mtx_file_path = args.s

    print(f'Success')
    print('Input:  ', input_file_path)
    print('Output: ', output_file_path)
    print('Matrix: ', mtx_file_path)


if __name__ == "__main__":
    main()
