import argparse
from pathlib import Path
import subprocess
from time import time

SAMP_SIZES = [1000000, 1500000, 2500000, 278197, 3000000, 5000000, 10000000, 30000000, 50000000, 100000000]
NUM_SAMPLES = [10, 5, 5, 5, 3, 3, 3, 3, 3, 3]


def parse_args():
    parser = argparse.ArgumentParser(
        prog="measure-performance",
        description="Measure performance of FASTX quality statistics executable"
    )
    parser.add_argument("-e", "--exec", help="Path to fastx quality statistics executable", type=str, required=True)
    parser.add_argument("-g", "--gen", help="Path to sample generator executable", type=str, required=True)
    parser.add_argument("-t", "--tmp", help="Temporary sample directory", type=str, required=True)

    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args):
    if not Path(args.exec).is_file():
        raise ValueError(f"Executable not found: {args.exec}")
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
    return elapsed_ms, result.stdout


def run_measurement(args):
    for i, sample_size in enumerate(SAMP_SIZES):
        size_str = str(sample_size)

        for session in range(NUM_SAMPLES[i]):
            generate_sample(args.gen, args.tmp, size_str)
            elapsed_ms, _ = execute_quality_stats(args.exec, args.tmp, size_str)
            print(f"[{size_str}:{session}] Time: {elapsed_ms:.2f}ms")
            delete_sample(args.tmp, size_str)


def main():
    try:
        args = parse_args()
        run_measurement(args)
    except Exception as e:
        print(f"Error: {e}")
        raise


if __name__ == "__main__":
    main()
