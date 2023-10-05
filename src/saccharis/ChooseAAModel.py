###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
import atexit
import glob
import os
import subprocess
from collections import defaultdict
from subprocess import PIPE
import random
import sys
from enum import Enum
from functools import reduce
from logging import getLogger, Logger

from Bio import SeqIO
from saccharis.CazyScrape import Mode
from saccharis.Muscle_Alignment import main as muscle
from saccharis.utils.PipelineErrors import AAModelError
from saccharis.utils.Formatting import convert_path_wsl


class TreeBuilder(Enum):
    RAXML = 1
    FASTTREE = 2
    RAXML_NG = 3


def compute_subsample(pruned_list, family, output_folder, num_threads, scrape_mode):
    subsample_size = 4000
    # Create Directory for muscle and change to this directory
    sub_folder = os.path.join(output_folder, "subsample")
    sub_file = os.path.join(sub_folder, f"{family}_subsample.fasta")
    if not os.path.isdir(sub_folder):
        os.mkdir(sub_folder)

    print("########################################################################################")
    print("########################################################################################")
    print("# Creating a Subsample of the pruned sequences for amino acid modeling as the original #")
    print("#      fasta file has over 4000 sequences, which is beyond the bounds of prottest      #")
    print("########################################################################################")
    print("########################################################################################\n")

    # Extract the Subsample of Sequences
    print("Extracting subsample...\n\n")

    # cazy_file = family + "_subsample.fasta"
    # $cmd1 = "fasta_subsample.pl $pfile 1500 -seed 7 > $cazy_file; "
    # &run_cmd($cmd1, $cmd2)
    # $cazy_file = $local . "/" . $cazy_file
    # Call dbcan
    # print("Running dbcan...\n\n")
    # dbcan_file = dbcan(cazy_file, family, local)
    # todo: QUESTION: why don't we just take 4000 dbcan pruned seqs and then run muscle on that? This is basically the
    #  approach I have used, since it saves time, no point recomputing hmmers.

    #   takes the list of seqs in BioSeq format and writes a random sample of 4000 to file
    random.seed("SACCHARIS", 2)  # Why yes, I AM using the program name and major revision number as a random seed
    sample_seqs = random.sample(pruned_list, subsample_size)
    with open(sub_file, 'w') as file:
        SeqIO.write(sample_seqs, file, 'fasta')

    # Call muscle
    # print("Running muscle...\n\n")
    paths = muscle(sub_file, subsample_size, family, scrape_mode, sub_folder, id_dict=None, is_subsample=True,
                   threads=num_threads)
    muscle_subsample_path = paths[2]

    print("#######################################################################################")
    print("#######################################################################################")
    print("#  Subsample process has been created and returning new subsample muscle filename     #")
    print("#   for prottest testing                                                              #")
    print("#######################################################################################")
    print("#######################################################################################\n")

    return muscle_subsample_path


def parse_best_model(outpath, tree_program, modeltest=False):

    # Parse the prottest results to obtain the model for raxML
    raxmodel = defaultdict(int)
    # i, g = (False, False) # probably unnecessary, need to test
    models = []
    with open(outpath, 'r') as protfile:
        for line in protfile:
            if modeltest:
                # if line.__contains__("Best model according to"):
                if line[0:6] == "Model:":
                    assert (line[20:-1] == line.split(' ')[-1].strip())  # ensure new method equals old, delete later
                    models.append(line.split(' ')[-1].strip())
            else:
                if line[0:4] == "Best":
                    models.append(line[line.find(':')+2:-1])

    # Use models parsed from the file to create the raxml modelname and push to dict incrementing identical values
    for model in models:
        # todo: this section is vastly simplified for fasttree, but can we capture number after the G? does fasttree
        #  support that? What about +R? Answer: fasttree does not, but maybe old raxml 8 does?
        # todo: move the parsing of model type for fasttree and raxml into respective tree build functions. Would
        #  simplify adding BLOSUM support to fasttree by packaging a BLOSUM matrix file with SACCHARIS, and makes
        #  extending to other tree software easier in the future.
        i = model.__contains__("+I")
        g = model.__contains__("+G")
        matrix = model.split('+')[0]

        # Set Tree ModelName based on RAxML or FastTree
        if tree_program == TreeBuilder.RAXML:   # Create the RaxML ModelName
            rxm = "PROT" + ("GAMMA" if g else "CAT") + ("I" if i else "") + matrix
        elif tree_program == TreeBuilder.FASTTREE:   # Create the FastTree ModelName
            rxm = ("gamma" if g else "cat") + \
                  ("-wag" if matrix == "WAG" else "-lg" if matrix == "LG" else "-jtt")
        elif tree_program == TreeBuilder.RAXML_NG:
            rxm = model
        else:
            raise AAModelError(f"Invalid tree program specified: {tree_program} is not valid. This is some kind of "
                               f"internal bug, please report it to the developer through github or email.")

        # i, g = (False, False) # probably unnecessary, need to test

        # increment model count, as identified by hash
        raxmodel[rxm] += 1

    # Set - final matrix name to the hash key with the largest count
    best_tree_model = reduce(lambda a, b: a if raxmodel[a] > raxmodel[b] else b, raxmodel)

    return best_tree_model


def compute_best_model(muscle_input_file, pruned_list, family, output_folder, number_seqs, scrape_mode, num_threads=4,
                       tree_program=TreeBuilder.FASTTREE, force_update=False, user_run=None,
                       prottest_folder="/usr/local/prottest-3.4.2", use_modelTest=True, logger: Logger = getLogger()):
    # Create Directory for prottest and change to this directory
    if user_run:
        prot_file_path = os.path.join(output_folder, f"{family}_{tree_program.name}_{user_run:05d}model.txt")
    else:
        prot_file_path = os.path.join(output_folder, f"{family}_{tree_program.name}model.txt")

    # Need to confirm file has not been made already and if so get model from file
    if os.path.isfile(prot_file_path) and not force_update:
        print("Mutation modeling step has been run already, continuing with script...\n")
        with open(prot_file_path, 'r') as f:
            best_tree_model = f.readline()
    elif number_seqs < 3:  # prottest triggers an exception with less than 3 sequences
        # Use gamma distribution with wag matrix as default, since this is a common best model result
        if tree_program == TreeBuilder.RAXML:   # Default RaxML ModelName
            best_tree_model = "PROTGAMMAWAG"
        else:   # Default FastTree ModelName
            best_tree_model = "gamma-wag"

        with open(prot_file_path, 'w') as file:
            file.write(best_tree_model)
        print(f"INFO: Could not run mutation modelling with less than 3 sequences, so assuming {best_tree_model} model")
    else:
        if force_update:
            # check for partial run files and delete them here
            files_to_delete = glob.glob(os.path.join(output_folder, '*'))
            for file_path in files_to_delete:
                try:
                    os.remove(file_path)
                    logger.info(f"Deleted {file_path}")
                except OSError as e:
                    logger.info(f"Error deleting {file_path}: {e}")

        # Before Actually running prottest we need to make sure there are not more than 4000 sequences in the input
        # file, prottest will not run with more than 4000 sequences, we need to make a subsample dataset to run with
        subsample_file = os.path.join(output_folder, "subsample", family + ".muscle_aln.phyi")
        if os.path.isfile(subsample_file):
            print("Subsample step has been run prior, Continuing\n\n")
            muscle_input_file = subsample_file
        else:
            if number_seqs >= 4000:
                # pass in pruned list to avoid recomputing hmmer bounds
                muscle_input_file = compute_subsample(pruned_list, family, output_folder, num_threads, scrape_mode)

        # todo: figure out why this section appears to be unnecessary
        # Copy the properties file for prottest to the cwd
        # my $prottest_properties = output_folder . "/prottest.properties";
        # if ( -f $prottest_properties ):
        #     print "Properties File has already been copied, proceeding with prottest run\n\n";
        # else:
        #     $cmd1 = "cp /usr/local/prottest3/prottest.properties .; ";
        #     &run_cmd($cmd1, $cmd2);
        #     $cmd1 = "";

        # Run prottest on the muscle results
        print(f"Running {'modeltest-ng' if use_modelTest else 'prottest'} - search for best model\n")
        modeltest_outfile = family + ".modeltest.out"
        modeltest_outpath = os.path.join(output_folder, modeltest_outfile)
        outfile = family + ".prottest.out"
        outpath = os.path.join(output_folder, outfile)

        # todo: load prottest_path from environment variable? Since it's not on $PATH, user may need to specify
        #  custom location. If this is changed, add a check that the file exists else raise a prompt to get location?
        #  maybe always download prottest to install folder and use it from there?
        prottest_path = os.path.join(prottest_folder, "prottest-3.4.2.jar")
        if not os.path.isfile(prottest_path) and not use_modelTest:
            print(f"prottest-3.4.2.jar not found at path: {prottest_path}\n"
                  f"Please install it! There is a linux install script in the SACCHARIS 2 install directory!\n"
                  f"Exiting...")
            exit()

        if use_modelTest:
            #     todo: pick starting topology for modeltest?
            if sys.platform.startswith("win"):
                try:
                    win_muscle_input_file = convert_path_wsl(muscle_input_file)
                    win_outpath = convert_path_wsl(modeltest_outpath[0:-4])

                    args = ["wsl", "modeltest-ng", "-d", "aa", "-i", win_muscle_input_file, "-o", win_outpath, "-h",
                            "uigf", "-p", f"{num_threads}"]
                except subprocess.CalledProcessError as error:
                    print(error.args)
                    print("wsl and/or wslpath are not installed. Please install or update windows subsystem for linux"
                          "to the latest version.\nExiting SACCHARIS...")
                    sys.exit(2)

            else:
                args = ["modeltest-ng", "-d", "aa", "-i", muscle_input_file, "-o", modeltest_outpath[0:-4], "-h",
                        "uigf", "-p", f"{num_threads}"]

            if tree_program == TreeBuilder.FASTTREE:
                args += ["-m", "JTT,LG,WAG"]

            if force_update:
                args += ["--force"]
        else:
            # todo:remove old prottest code
            args = ["java", "-jar", prottest_path, "-i", muscle_input_file, "-o", outpath, "-S", "0",
                    "-all-distributions", "-AIC", "-AICC", "-BIC", "-DT", "-threads", f"{num_threads}"]
            if tree_program == TreeBuilder.FASTTREE:
                args += ["-JTT", "-LG", "-WAG"]

        try:
            # subprocess.run(args, check=True)
            msg = f"modeltest args: {' '.join(args)}"
            logger.debug(msg)
            # main_proc = subprocess.Popen(args, text=True, stderr=PIPE, stdout=PIPE)
            main_proc = subprocess.Popen(args, text=True, stderr=PIPE)
            atexit.register(main_proc.kill)
            outstring, errstring = main_proc.communicate()
            # main_proc.wait()
            atexit.unregister(main_proc.kill)
            if main_proc.returncode != 0:
                retcode_msg = f"modeltest process did not return correctly. Process return code: {main_proc.returncode}"
                # logger.debug(f"modeltest STDOUT: {outstring}")
                logger.error(f"modeltest STDERR: {errstring}")
                logger.error(retcode_msg)
                raise AAModelError(retcode_msg)
        except subprocess.CalledProcessError as error:
            raise AAModelError("modeltest failed for an unknown reason. Check that it is installed. If you are using"
                               "a windows Operating system, note that modeltest-ng is not compatible with windows, and"
                               "SACCHARIS tries to run it on the default windows subsystem for linux distribution and "
                               "you should install WSL, then install modeltest-ng on your WSL install.") from error

        print("modeltest finished - proceeding with parsing modeltest results\n")

        if use_modelTest:
            best_tree_model = parse_best_model(modeltest_outpath, tree_program, modeltest=use_modelTest)
        else:
            best_tree_model = parse_best_model(outpath, tree_program, modeltest=use_modelTest)

        tp_name = "FastTree" if tree_program == TreeBuilder.FASTTREE \
            else "RAxML-NG" if tree_program == TreeBuilder.RAXML_NG \
            else "RAxML"
        print(f"Parsing of output files completed and have best model: {best_tree_model} for {tp_name} Run\n\n")

        with open(prot_file_path, 'w') as file:
            file.write(best_tree_model)

    return best_tree_model


if __name__ == "__main__":
    test_family = "PL9"
    test_mode = Mode.CHARACTERIZED

    group_folder = os.path.join(os.getcwd(), "../../output", test_family, test_mode.name)
    input_file_path = os.path.join(group_folder, "muscle", "PL9_CHARACTERIZED.muscle_aln.phyi")
    test_output_folder = os.path.join(group_folder, "prottest")
    if not os.path.isdir(test_output_folder):
        os.mkdir(test_output_folder, 0o755)

    compute_best_model(input_file_path, "dbcan", "PL9", test_mode, test_output_folder, 13, num_threads=12)
