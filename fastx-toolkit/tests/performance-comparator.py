import argparse
import os
from pathlib import Path
import subprocess
from time import time

SAMP_SIZES = [1000000, 1500000, 2500000, 278197, 3000000, 5000000, 10000000, 30000000, 50000000, 100000000]
NUM_SAMPLES = [10, 5, 5, 5, 3, 3, 3, 3, 3, 3]


def parse_args():
    parser = argparse.ArgumentParser(
        prog="performance-comparator",
        description="Compare performance between old and new FASTX quality statistics implementations"
    )
    parser.add_argument("-o", "--old", help="Path to old fastx quality statistics executable", type=str, required=True)
    parser.add_argument("-n", "--new", help="Path to new fastx quality statistics executable", type=str, required=True)
    parser.add_argument("-g", "--gen", help="Path to sample generator executable", type=str, required=True)
    parser.add_argument("-t", "--tmp", help="Temporary sample directory", type=str, required=True)

    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args):
    if not Path(args.old).is_file():
        raise ValueError(f"Old executable not found: {args.old}")
    if not Path(args.new).is_file():
        raise ValueError(f"New executable not found: {args.new}")
    if not Path(args.gen).is_file():
        raise ValueError(f"Generator executable not found: {args.gen}")
    if not Path(args.tmp).is_dir():
        raise ValueError(f"Temporary directory not found: {args.tmp}")


def generate_sample(gen_path, tmp_dir, size):
    sample_path = Path(tmp_dir) / str(size)
    command = f'{gen_path} -s fastq --nr {size} -o {sample_path} --mns 50 --mxs 200'
    subprocess.run(command, capture_output=False, text=False, check=True, shell=True)


def delete_sample(tmp_dir, size):
    sample_path = Path(tmp_dir) / str(size)
    sample_path.unlink(missing_ok=True)


def execute_quality_stats(program_path, tmp_dir, size):
    sample_path = Path(tmp_dir) / str(size)
    command = f'{program_path} -i {sample_path}'

    start_time = time()
    result = subprocess.run(command, capture_output=True, text=True, check=True, shell=True)
    end_time = time()

    elapsed_ms = (end_time - start_time) * 1000
    output = result.stdout.replace('\n', '').replace('\r', '').strip()

    return elapsed_ms, output


def compare_results(old_result, new_result, size, session):
    old_time, old_output = old_result
    new_time, new_output = new_result

    status = 'PASSED' if old_output == new_output else 'MISMATCH'
    print(f"[{size}:{session}] {status}, Old: {old_time:.2f}ms, New: {new_time:.2f}ms")


def run_comparison(args):
    for i, sample_size in enumerate(SAMP_SIZES):
        size_str = str(sample_size)
        generate_sample(args.gen, args.tmp, size_str)

        for session in range(NUM_SAMPLES[i]):
            old_result = execute_quality_stats(args.old, args.tmp, size_str)
            new_result = execute_quality_stats(args.new, args.tmp, size_str)
            compare_results(old_result, new_result, size_str, session)

        delete_sample(args.tmp, size_str)


def main():
    try:
        args = parse_args()
        run_comparison(args)
    except Exception as e:
        print(f"Error: {e}")
        raise


if __name__ == "__main__":
    main()
