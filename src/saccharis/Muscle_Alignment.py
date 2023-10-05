###############################################################################
# Perform multiple sequence alignment for tree building using muscle v3 or v5
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
# Built in libraries
import csv
import math
import os
import re
import subprocess
import json
import sys
import atexit

# Dependency libraries
import Bio.SeqIO

# Internal imports
from saccharis.CazyScrape import Mode
from saccharis.utils.PipelineErrors import PipelineException


def muscle_rename(infile, outfile, id_dict):
    renamed_muscle_data = ""
    with open(infile, 'r', newline='\n') as infile:
        muscle_data = csv.reader(infile, delimiter=' ')
        for line in muscle_data:
            if not line:  # add newline for empty line
                renamed_muscle_data += '\n'
            else:
                if line[0] in id_dict:
                    # replace 10 character temporary id with original genbank accession from dict
                    # note that this does NOT replace the user sequence IDs of the form "U*********" in SACCHARIS,
                    # because those IDs are not put in this id_dict file
                    line[0] = id_dict[line[0]]
                renamed_muscle_data += ' '.join(line) + '\n'

    with open(outfile, 'w') as ren_file:
        ren_file.write(renamed_muscle_data)

    return outfile


def main(input_file, seq_count, family, mode, output_folder, id_dict, force_update=False, is_subsample=False,
         user_run_id=None, threads=math.ceil(os.cpu_count() * 0.75)):

    # setup output folder and file names
    if user_run_id is not None:
        muscle_path = os.path.join(output_folder, f"{family}_{mode.name}_UserRun{user_run_id:05d}.muscle_aln.phyi")
    else:
        muscle_path = os.path.join(output_folder, f"{family}_{mode.name}.muscle_aln.phyi")
    muscle_path_efa = re.sub(r"aln\.phyi", "aln.efa", muscle_path)
    muscle_ren_path = re.sub(r"aln\.phyi", "aln_mod.phyi", muscle_path)
    muscle_fast = re.sub(r"aln\.phyi", "aln_mod_fast.phyi", muscle_path)
    if not os.path.exists(output_folder):
        os.mkdir(output_folder, 0o755)

    if os.path.isfile(muscle_path_efa) and os.path.getsize(muscle_path_efa) == 0:
        # muscle writes an empty file, even if a prior run failed or was interrupted
        os.remove(muscle_path_efa)

    # Run Muscle Alignment
    if os.path.exists(muscle_path) and force_update is False:
        print("Muscle has already been run, continuing script\n")
    elif seq_count < 2:
        # We need to convert the file instead of running muscle because muscle doesn't write the output file when it
        # doesn't perform an alignment with a single sequence
        print("Less than two input sequences to muscle, no reason to align!")
        print("Skipping muscle alignment...")
        Bio.SeqIO.convert(input_file, 'fasta', muscle_path, 'phylip')
    else:
        try:
            print("Running the muscle alignment on the pruned fasta data\n")
            # get muscle version
            proc_out = subprocess.run(["muscle", "-version"], capture_output=True, text=True, check=True)
            if proc_out.stdout.__contains__("MUSCLE v3"):
                args = ["muscle", "-in", input_file, "-phyiout", muscle_path]
            elif proc_out.stdout.__contains__("muscle v5") or proc_out.stdout.__contains__("muscle 5"):
                args = ["muscle", "-align", input_file, "-output", muscle_path_efa, "-threads", str(threads)]
                # subprocess.run(args, check=True)
            else:
                print("Unrecognized version of muscle, aborting...")
                sys.exit(1)
            # subprocess.run(args, check=True)
            main_proc = subprocess.Popen(args)
            atexit.register(main_proc.kill)
            main_proc.wait()
            atexit.unregister(main_proc.kill)
            if main_proc.returncode != 0:
                raise PipelineException("Muscle alignment process failed to return properly.")
            if proc_out.stdout.__contains__("muscle v5") or proc_out.stdout.__contains__("muscle 5"):
                Bio.SeqIO.convert(muscle_path_efa, "fasta", muscle_path, "phylip")
            print("Muscle Alignment completed\n\n")
        except subprocess.CalledProcessError as error:
            raise PipelineException("Muscle could not execute properly.") from error
        finally:
            if os.path.isfile(muscle_path_efa) and os.path.getsize(muscle_path_efa) == 0:
                # muscle writes an empty file even when it fails
                os.remove(muscle_path_efa)

    # Only need muscle_path for subsampling to choose aa model
    if is_subsample:
        return None, muscle_path, None

    # Run the Muscle Renamer script
    if os.path.exists(muscle_ren_path) and force_update is False:
        print("Renamed muscle file already exists, continuing script")
    else:
        print("Renaming muscle output...")
        muscle_rename(muscle_path, muscle_ren_path, id_dict)
        print("Muscle Rename Completed!")

    # Create the FastTree Variant of muscle file
    if os.path.exists(muscle_fast) and force_update is False:
        print("Muscle - FastTree Variant has been created, continuing script")
    else:
        print("Creating a FastTree variant file of the muscle output")
        # Original perl parsing line
        # $cmd1 = "awk '{if (\$1 \~ \/[0-9]\/) print \$0 \ else print \"\",\$0}' $muscle_ren_path | sed 's/^ \$//' >
        # $muscle_fast " &run_cmd($cmd1, $cmd2)

        # Python subprocess call for unix systems
        # subprocess.run("awk '{if ($1 ~ /[0-9]/) print $0; else print \"\",$0}' " + muscle_ren_path +
        #                " | sed 's/^ $//' > " + muscle_fast, shell=True, check=True)

        # Above approaches have been replaced by a python native conversion.
        # It seems to just adds spaces to every line that isn't a sequence start, which is detemined by containing a
        # number in the first space delimited string. I assume that this is to line up the character positions of amino
        # acid sequences, presumably because fasttree indexes the subsequent amino acid characters by skipping the same
        # number of characters seen in the first block.

        try:
            num_regex = re.compile("\d+")
            text = ""
            with open(muscle_ren_path, 'r') as fast_in_file:
                for line in fast_in_file:
                    try:
                        if line == '\n' or num_regex.search(line.split(' ')[0]) or num_regex.search(line.split(' ')[1]):
                            text += line
                        else:
                            text += ' ' + line
                    except IndexError:
                        text += ' ' + line

            with open(muscle_fast, 'w') as fast_out_file:
                fast_out_file.write(text)
        except IOError as error:
            raise PipelineException("Could not parse muscle output properly.") from error

        print("FastTree Variant Created")

    # The plain muscle file is easiest for prottest so will us it with subsample run
    # if run == 0:
    #     return muscle_ren_path
    # else:
    #     return muscle_path
    # Instead of above, we return all 3 file paths and let the calling function decide which one it wants instead of
    # using the run argument
    return muscle_ren_path, muscle_path, muscle_fast


if __name__ == "__main__":
    test_family = "PL9"
    test_mode = Mode.CHARACTERIZED

    group_folder = os.path.join(os.getcwd(), "../../output", test_family, test_mode.name)
    input_file_path = os.path.join(group_folder, "dbcan2", f"{test_family}_CHARACTERIZED_cazy.mod.pruned.fasta")
    test_output_folder = os.path.join(group_folder, "muscle")
    id_file = os.path.join(group_folder, "dbcan2", f"{test_family}_{test_mode.name}_key_id_pairs.json")
    with open(id_file, 'r', encoding='utf-8') as f:
        test_id_dict = json.loads(f.read())
    main(input_file_path, 0, "PL9", Mode.CHARACTERIZED, test_output_folder, test_id_dict)
