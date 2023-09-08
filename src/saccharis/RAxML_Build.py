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
from logging import Logger, getLogger
from math import ceil
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
    # muscle_depth = `head -1 $muscle` # muscle depth = first line in muscle file
    # @md = split(/ /, $muscle_depth)
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
                        user_run: int = None, logger: Logger = getLogger()):
    if user_run:
        # todo: probably need to update these filenames
        rax_tree = f"RAxML_bipartitions.{user_run:05d}.A1"
    else:
        rax_tree = "RAxML_bipartitions.A1"
    bootstrap = 100
    file_output_path = os.path.join(output_dir, rax_tree)

    return file_output_path


if __name__ == "__main__":
    test_family = "PL9"
    test_mode = "CHARACTERIZED"
    group = os.path.join(os.getcwd(), "../../output", test_family, test_mode)
    muscle_input = os.path.join(group, "muscle", "PL9_CHARACTERIZED.muscle_aln_mod.phyi")
    raxml_folder = os.path.join(group, "raxml")
    if not os.path.isdir(raxml_folder):
        os.mkdir(raxml_folder, 0o755)
    main(muscle_input, "PROTGAMMAIWAG", raxml_folder, "raxmlHPC-PTHREADS-AVX", 13, threads=8)
