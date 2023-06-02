###############################################################################
# FastTree - Use muscle alignment and the model chosen by modeltest to build
#             tree - faster than RAxML
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
# Built in libraries
import os
import subprocess
import atexit

# Internal Imports
from saccharis.utils.PipelineErrors import PipelineException


def main(muscle_input_path, amino_model, output_dir, force_update=False, user_run=None):
    out_filename = f"FastTree_bootstrap_UserRun{user_run:05d}.tree" if user_run else f"FastTree_bootstrap.tree"
    output_file_path = os.path.join(output_dir, out_filename)
    # parse the amino acid model to set proper flags in run
    model = amino_model.split('-')

    # Set args and run fasttree
    if os.path.isfile(output_file_path) and not force_update:
        print("\n\nFastTree has built this tree before, loading tree data from previous run.\n\n")
        return output_file_path
    else:
        print("Building best tree - using FastTree\n")
        model[0] = "" if model[0] == "cat" else "-gamma "
        model[1] = "-wag " if model[1] == "wag" else "-lg " if model[1] == "lg" else ""
        # todo: set option to allow use of the multithreaded but non-determistic version of fasttree?
        # Some repos use "fasttree" as command name, others use "FastTree", so we test what's installed
        try:
            subprocess.run("fasttree -help", shell=True, capture_output=True, check=True)
            command_name = "fasttree"
        except subprocess.CalledProcessError:
            try:
                subprocess.run("FastTree -help", shell=True, capture_output=True, check=True)
                command_name = "FastTree"
            except subprocess.CalledProcessError:
                raise UserWarning("fasttree is not installed! Make sure it is available on path via the 'fasttree' "
                                  "command")

        command = f"{command_name} {model[1]}{model[0]}-out {output_file_path} {muscle_input_path}"

        # proc_out = subprocess.run(command, shell=True)
        proc_out = subprocess.Popen(command, shell=True)
        atexit.register(proc_out.kill)
        proc_out.wait()
        atexit.unregister(proc_out.kill)

        if proc_out.returncode == 0:
            print("FastTree has finished\n\n")
        else:
            raise PipelineException("ERROR: FastTree did not return valid output. Exiting pipeline...")

    return output_file_path


if __name__ == "__main__":
    test_family = "PL9"
    test_mode = "CHARACTERIZED"
    group = os.path.join(os.getcwd(), "../../output", test_family, test_mode)
    fasttree_folder = os.path.join(group, "fasttree")
    if not os.path.isdir(fasttree_folder):
        os.mkdir(fasttree_folder, 0o755)
    muscle_input = os.path.join(group, "muscle", "PL9_Mode.CHARACTERIZED.muscle_aln_mod_fast.phyi")

    with open(os.path.join(group, "prottest", "PL9_prot_model.txt")) as f:
        aa_model = f.readline()

    main(muscle_input, aa_model, fasttree_folder)
