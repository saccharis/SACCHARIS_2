###############################################################################
# RaxML - Use muscle alignment and the model chosen by modeltest to build
#          the best tree. Slower than fasttree, possible higher quality.
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
# Built in libraries
import atexit
import glob
import os
import subprocess
import sys
from logging import Logger, getLogger
from math import ceil

from saccharis.utils.Formatting import convert_path_wsl
# Internal Imports
from saccharis.utils.PipelineErrors import PipelineException


def main(muscle_input_file: str | os.PathLike, amino_model: str, output_dir: str | os.PathLike,
         raxml_version: str, num_seqs: int, threads: int = 4, force_update: bool = False,
         user_run: int = None, logger: Logger = getLogger()):
    if user_run:
        rax_tree = f"RAxML_bipartitions.{user_run:05d}.A1"
    else:
        rax_tree = "RAxML_bipartitions.A1"
    bootstrap = 100

    # Calculate an optimal number of threads to use
    opt_thr = ceil(num_seqs / 300)  # optimal thread num = # seqs / 300 rounded up
    # todo: QUESTION: why do we do this? alignments of smaller groups still seem faster with more threads

    # Now set threads used to the lower value of $opt_thr and $threads
    #  todo: should we up the thread count above 16? Why is it limited to 16 threads anyway?
    if threads < opt_thr:
        opt_thr = threads
    elif threads >= 16:
        opt_thr = 16
    elif opt_thr == 1:
        opt_thr = 2  # raxML pthread version requires at least two threads

    file_output_path = os.path.join(output_dir, rax_tree)

    if os.path.isfile(file_output_path) and not force_update:
        msg = "\n\nRAxML has built this tree before, loading tree data from previous run.\n\n"
        print(msg)
        logger.debug(msg)
        return file_output_path
    else:

        # check for partial run files and delete them here
        files_to_delete = glob.glob(os.path.join(output_dir, '*.A1'))
        for file_path in files_to_delete:
            try:
                os.remove(file_path)
                logger.info(f"Deleted {file_path}")
            except OSError as e:
                logger.info(f"Error deleting {file_path}: {e}")

        print("Building best tree - using RAxML\n")
        # thread_arg = "" if opt_thr == 1 else f"-T {opt_thr} "
        extension = f"UserRun{user_run:05d}.A1" if user_run else "A1"
        command = f"{raxml_version} -f a -T {opt_thr} -p 0111 -s {muscle_input_file} -n {extension} " \
                  f"-m {amino_model} -x 0123 -# {bootstrap}"

# if ($opt_thr == 1) {
        #     $cmd1 = f"{raxml_version} -f a -p 0111 -s $muscle -n A1 -m $tree -x 0123 -# $bootstrap; ";
        # } else {
        #     $cmd1 = "$raxml -f a -T $opt_thr -p 0111 -s $muscle -n A1 -m $tree -x 0123 -# $bootstrap; ";
        # }
        # &run_cmd($cmd1, $cmd2);

        msg = f"Running command: {command}"
        logger.info(msg)
        # subprocess.run(command, shell=True, cwd=output_dir, check=True)
        main_proc = subprocess.Popen(command, shell=True, cwd=output_dir)
        atexit.register(main_proc.kill)
        main_proc.wait()
        atexit.unregister(main_proc.kill)
        if main_proc.returncode != 0:
            msg = f"raxml tree building process failed to return properly. Return code: {main_proc.returncode}"
            logger.error(msg)
            raise PipelineException(msg)

        msg = "RaxML has finished\n\n"
        print(msg)
        logger.debug(msg)
        return file_output_path


def build_tree_raxml_ng(muscle_input_file: str | os.PathLike, amino_model: str, output_dir: str | os.PathLike,
                        num_seqs: int, threads: int = 4, force_update: bool = False,
                        user_run: int = None, logger: Logger = getLogger(), bootstraps: int = 100):
    if user_run:
        # todo: update these prefixes to be more informative? Is this redudant?
        rax_tree = f"RAxML-NG_output.U{user_run:05d}"
    else:
        rax_tree = "RAxML-NG_output"

    initial_seed = 111
    file_output_prefix = os.path.join(output_dir, rax_tree)

    if sys.platform.startswith("win"):
        command = ["wsl"]
        muscle_input_path = convert_path_wsl(muscle_input_file)
        file_output_path = convert_path_wsl(file_output_prefix)
        validity_args = ["wsl"]
    else:
        command = []
        muscle_input_path = muscle_input_file
        file_output_path = file_output_prefix
        validity_args = []

    validity_args += ["raxml-ng", "--parse", "--seed", str(initial_seed), "--msa", muscle_input_path, "--model", amino_model, "--prefix", file_output_path]
    try:
        valid_result = subprocess.run(validity_args, capture_output=True, encoding="utf-8", check=True)
        logger.info(valid_result.stdout)
    except FileNotFoundError as err:
        logger.error(err)
        msg = "raxml-ng not found, check that it's available via PATH variable."
        logger.error(msg)
        raise PipelineException(msg) from err

    optimal_threads = threads
    can_run = False
    for line in valid_result.stdout.split('\n'):
        if line.__contains__("Recommended number of threads"):
            optimal_threads = int(line.split(' ')[-1])
        elif line.__contains__("Alignment can be successfully read"):
            can_run = True

    if not can_run:
        logger.error(valid_result.stdout)
        logger.error("RAxML-NG cannot read MSA.")
        raise PipelineException("RAxML-NG cannot read MSA.")
    run_threads = min(optimal_threads, threads)

    # todo: add outgroup option --outgroup [csv list]
    command += ["raxml-ng", "--all", "--threads", f"auto{'{' + str(run_threads) + '}'}", "--seed", str(initial_seed),
                "--msa", muscle_input_path, "--prefix", file_output_path, "--model", amino_model,
                "--bs-trees", str(bootstraps)]

    if force_update:
        command += ["--redo"]

    msg = f"Running command: {command}"
    logger.info(msg)

    try:
        assert(os.path.isfile(muscle_input_file))
        assert(os.path.isdir(output_dir))
        main_proc = subprocess.Popen(command, cwd=output_dir)
        atexit.register(main_proc.kill)
        main_proc.wait()
        atexit.unregister(main_proc.kill)
        if main_proc.returncode != 0:
            msg = f"raxml-ng tree building process failed to return properly. Return code: {main_proc.returncode}"
            logger.error(msg)
            raise PipelineException(msg)

    except FileNotFoundError as err:
        logger.error(err)
        msg = "raxml-ng not found, check that it's available via PATH variable."
        logger.error(msg)
        raise PipelineException(msg) from err

    # todo: check if bootstraps converged and then run until they converge. This may take large amounts of compute and
    #  should be an optional feature. Will have to recompute best tree with support, so maybe just do this
    #  earlier instead of running with --all?
    # see https://github.com/amkozlov/raxml-ng/wiki/Tutorial#bootstrapping
    # bsconverge_args = ["raxml-ng", "--bsconverge", "--bs-trees", f"{file_output_path}.raxml.bootstraps",
    #                    "--prefix", f"{file_output_path}", "--seed", "2", "--threads", "2", "--bs-cutoff", "0.01"]
    # if sys.platform.startswith("win"):
    #     bsconverge_args.insert(0, "wsl")
    # bsconverge_result = subprocess.run(bsconverge_args, capture_output=True, encoding="utf-8", check=True)

    msg = "RaxML-NG has finished\n\n"
    print(msg)
    logger.debug(msg)

    best_tree_with_support_path = os.path.join(f"{file_output_prefix}.raxml.support")
    return best_tree_with_support_path


if __name__ == "__main__":
    test_family = "PL9"
    test_mode = "CHARACTERIZED"
    group = os.path.join(os.getcwd(), "../../output", test_family, test_mode)
    muscle_input = os.path.join(group, "muscle", "PL9_CHARACTERIZED.muscle_aln_mod.phyi")
    raxml_folder = os.path.join(group, "raxml")
    if not os.path.isdir(raxml_folder):
        os.mkdir(raxml_folder, 0o755)
    main(muscle_input, "PROTGAMMAIWAG", raxml_folder, "raxmlHPC-PTHREADS-AVX", 13, threads=8)
