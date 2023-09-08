###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
# Built in libraries
import argparse
import csv
import math
import os
import sys
import re
import json
from csv import DictReader

# Dependency libraries
from dbcan_cli import run_dbcan
from Bio import SeqIO

# Internal imports
from saccharis.utils.PipelineErrors import UserError
from saccharis.utils.FamilyCategories import Matcher
from saccharis.utils.PipelineErrors import PipelineException
from saccharis.utils.AdvancedConfig import get_db_folder
from saccharis.utils.DatabaseDownload import download_database


def filter_prune(fasta_filepath, bounds_file, family, output_folder, source, prune=True, accession_dict=None):
    #   Filter hmmer output for correct family and unique accession numbers
    with open(bounds_file, 'r', newline='\n') as hmmer_tsv:
        if source == "dbcan":
            entry_reader = csv.reader(hmmer_tsv, delimiter='\t')
            hmmer_list = list(entry_reader)
            family_column = 0
            accession_column = 2
            gene_start_column = 7
            gene_end_column = 8
            matcher = Matcher()
            hmmer_list = [entry for entry in hmmer_list if entry[0] != "HMM Profile"]
            for entry in hmmer_list:
                # family_extracted = matcher.extract_cazy_family(entry[family_column])
                # if family_extracted == "" or None:
                #     raise PipelineException(f"Bad family info for entry in dbcan hmmer output. Element info:{entry}")
                # entry[family_column] = family_extracted
                entry[family_column] = entry[family_column][0:entry[family_column].find('.')]
        elif source == "pfam":
            lines = hmmer_tsv.readlines()
            lines = [line.strip() for line in lines if not line.__contains__('#') and len(line.strip()) > 1]
            lines = [re.sub(" +", " ", line) for line in lines]
            lines = [line.split(' ') for line in lines]
            hmmer_list = lines
            family_column = 5
            accession_column = 0
            gene_start_column = 1
            gene_end_column = 2
        else:
            raise UserError("Wrong output source of data file to filter")

    # filters correct family
    if family.__contains__('_') or source == "pfam":
        # match full family and subfamily
        hmmer_list_filtered = [entry for entry in hmmer_list if entry[family_column] == family]
    else:
        # match family only
        hmmer_list_filtered = [entry for entry in hmmer_list if entry[family_column].split('_')[0] == family]

    # initialize counts for each unique accession #
    hmmer_counts = {entry[accession_column]: 0 for entry in hmmer_list_filtered}
    hmmer_renamed_count = {}
    #   Count occurrences of each accession in hmmer output
    for entry in hmmer_list_filtered:
        hmmer_counts[entry[accession_column]] += 1
        if hmmer_counts[entry[accession_column]] > 1:
            hmmer_renamed_count[entry[accession_column]] = 0  # initialize state variable for renaming

    if prune:
        #   Append "[<Occurrence count>]" to non-unique accession numbers
        for entry in hmmer_list:
            if entry[accession_column] in hmmer_counts and hmmer_counts[entry[accession_column]] > 1:
                hmmer_renamed_count[entry[accession_column]] += 1
                entry[accession_column] += f"<{hmmer_renamed_count[entry[accession_column]]}>"
    else:
        #   If not pruning, need to remove duplicate accessions, since we would then have two identical sequences
        for entry in hmmer_list:
            if entry[accession_column] in hmmer_counts and hmmer_counts[entry[accession_column]] > 1:
                hmmer_renamed_count[entry[accession_column]] += 1
                if hmmer_renamed_count[entry[accession_column]] > 1:
                    entry[accession_column] += f"<r>"  # mark for removal
        hmmer_list_filtered = [entry for entry in hmmer_list_filtered if not entry[accession_column].__contains__("<r>")]

    if source == "dbcan":
        hmmer_filename = re.sub(r"hmmer\.out", "hmmer_unique.out", os.path.basename(bounds_file))
    elif source == "pfam":
        hmmer_filename = re.sub(r"\.pfam", "_unique.pfam", os.path.basename(bounds_file))
    else:
        raise UserError("Wrong output source of data file to filter")

    #   Write filtered and unique hmmer output back to a file
    # todo: Perhaps delete this write? Doesn't seem that useful and I don't think it's ever needed
    if output_folder is not None:
        hmmer_outfile = os.path.join(output_folder, hmmer_filename)
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        with open(hmmer_outfile, 'w', newline='\n') as hmmer_tsv:
            entry_writer = csv.writer(hmmer_tsv, delimiter='\t')
            if len(hmmer_list_filtered) > 0:
                entry_writer.writerows(hmmer_list_filtered)
            else:
                print("# WARNING: No sequences matched the specified family! File not written.")

    #   prune sequences
    if accession_dict is None:
        accession_dict = SeqIO.to_dict(SeqIO.parse(fasta_filepath, "fasta"))
    pruned = []
    bounds_dict = {}
    #   note: hmmer_list_filtered is a SHALLOW COPY of hmmer_list, so it also has the "[<count>]" appended to names of
    #       non-unique accessions, since the sublists are shared between hmmer_list and hmmer_list_filtered
    for entry in hmmer_list_filtered:
        # entry_family = entry[0] # don't actually need this
        accession_unique = entry[accession_column]
        gene_start = entry[gene_start_column]
        gene_end = entry[gene_end_column]
        accession_short = re.sub(r"<\d+>", "", accession_unique)
        # Prune seq
        if prune:
            pruned_test = accession_dict[accession_short].seq[int(gene_start) - 1:int(gene_end)]
            bounds_dict[accession_unique] = (gene_start, gene_end)
        else:
            pruned_test = accession_dict[accession_short].seq
            bounds_dict[accession_unique] = (1, len(accession_dict[accession_short].seq)+1)
        pruned_desc = re.sub(accession_short, accession_unique, accession_dict[accession_short].description)
        pruned.append(SeqIO.SeqRecord(pruned_test, id=accession_unique, name=accession_unique,
                                      description=pruned_desc))

    return pruned, bounds_dict


def parse_eCAMI_dict(file_path):
    ecami_data = {}
    with open(file_path, 'r', newline='\n') as file:
        iterator = DictReader(file, dialect="unix", delimiter='\t')
        for row in iterator:
            ecami_data[row["protein_name"]] = {"fam_name:group_number": row["fam_name:group_number"],
                                               "subfam_name_of_the_group:subfam_name_count":
                                               row["subfam_name_of_the_group:subfam_name_count"].split('|')}
    return ecami_data


def parse_diamond_dict(file_path):
    diamond_data = {}
    with open(file_path, 'r') as file:
        iterator = DictReader(file, dialect="unix", delimiter="\t")
        for row in iterator:
            diamond_data[row["Gene ID"]] = {"CAZyme Predictions": row["CAZy ID"].split('|')}
    return diamond_data


def main(fasta_filepath, family, output_folder, mode, force_update=False, prune=True, accession_dict=None,
         threads=math.ceil(os.cpu_count() * .75), hmm_eval: float = 1e-15, hmm_cov: float = 0.35):

    download_database()

    # set up dbcan output filenames
    fasta_filename = os.path.split(fasta_filepath)[1]
    pruned_filepath = os.path.join(output_folder, re.sub(r"\.fasta", ".pruned.fasta", fasta_filename))
    fasta_mod_file = re.sub("pruned", "mod.pruned", pruned_filepath)
    id_file = os.path.join(output_folder, f"{family}_{mode.name}_key_id_pairs.json")
    bounds_file = os.path.join(output_folder, f"{family}_{mode.name}_bounds.json")
    ecami_file = os.path.join(output_folder, f"{family}_{mode.name}_eCAMI.json")
    diamond_file = os.path.join(output_folder, f"{family}_{mode.name}_diamond.json")

    try:
        if os.path.isfile(pruned_filepath) and os.path.isfile(fasta_mod_file) and os.path.isfile(id_file) \
                and os.path.isfile(bounds_file) and os.path.isfile(ecami_file) and os.path.isfile(diamond_file) \
                and not force_update:
            print("CAZymes already extracted, loading data from previous run.")
            print("If you would like to recalculate HMMERs, run SACCHARIS with --fresh")
            #  load and return existing data
            pruned = list(SeqIO.parse(pruned_filepath, 'fasta'))

            with open(id_file, 'r', encoding='utf-8') as f:
                mod_dict = json.loads(f.read())
            with open(bounds_file, 'r', encoding="utf-8") as f:
                bounds_dict = json.loads(f.read())
            with open(ecami_file, 'r', encoding="utf-8") as f:
                ecami_dict = json.loads(f.read())
            with open(diamond_file, 'r', encoding="utf-8") as f:
                diamond_dict = json.loads(f.read())

            return pruned, fasta_mod_file, mod_dict, bounds_dict, ecami_dict, diamond_dict
    except IOError as error:
        # todo: log error here
        print("Error loading data from previous run, recalculating instead...")

    # call dbcan2 script, which outputs results of HMMER, DIAMOND, and eCAMI to output_folder, if needed
    hmmer_filepath = os.path.join(output_folder, "hmmer.out")
    ecami_filepath = os.path.join(output_folder, "eCAMI.out")
    diamond_filepath = os.path.join(output_folder, "diamond.out")
    if not (os.path.exists(hmmer_filepath) and os.path.exists(ecami_filepath) and os.path.exists(diamond_filepath)) \
            or force_update:
        run_dbcan.run(fasta_filepath, "protein", outDir=output_folder, dbDir=get_db_folder(), dia_cpu=threads,
                      hmm_cpu=threads, tf_cpu=threads, stp_cpu=threads, eCAMI_jobs=threads, hmm_cov=hmm_cov,
                      hmm_eval=hmm_eval)

    print("Extract Sequences that hit family:", family)

    try:
        pruned, bounds_dict = filter_prune(fasta_filepath, hmmer_filepath, family, output_folder, "dbcan", prune,
                                           accession_dict)
        ecami_dict = parse_eCAMI_dict(ecami_filepath)
        diamond_dict = parse_diamond_dict(diamond_filepath)
    except KeyError as error:
        raise PipelineException("HMMER, Diamond, or eCAMI output files not in correct format. Have they changed their "
                                "csv output format? Please report this bug to the developer/maintainer through github!"
                                f"\nERROR MESSAGE: {error.args[0]}") from error

    # write pruned seqs to file
    with open(pruned_filepath, 'w', newline='\n') as f:
        SeqIO.write(pruned, f, 'fasta')

    # add modified id sequence id and write to file
    # line_list = fasta_data.split('\n')
    mod_dict = {}
    # mod_data = ""
    modified_count = 0
    for entry in pruned:
        # if entry.__contains__('>'):
        old_id = entry.id
        new_id = f"{modified_count:09d}"
        mod_dict[new_id] = old_id
        # mod_data += f">{modified_count:09d} {entry[1:-1]}\n"
        entry.id = new_id
        entry.name = new_id
        entry.description = new_id + " " + entry.description
        modified_count += 1
        # elif not entry == '':
        #     mod_data += entry + '\n'
    with open(fasta_mod_file, 'w', newline='\n') as f:
        # f.write(mod_data)
        SeqIO.write(pruned, f, 'fasta')

    # write dicts to translate modified ids back to genbank accessions to file and pass metadata to main pipeline
    with open(id_file, 'w', encoding='utf-8') as f:
        json.dump(mod_dict, f, ensure_ascii=False, indent=4)

    with open(bounds_file, 'w', encoding="utf-8") as f:
        json.dump(bounds_dict, f, ensure_ascii=False, indent=4)

    with open(ecami_file, 'w', encoding="utf-8") as f:
        json.dump(ecami_dict, f, ensure_ascii=False, indent=4)

    with open(diamond_file, 'w', encoding="utf-8") as f:
        json.dump(diamond_dict, f, ensure_ascii=False, indent=4)

    return pruned, fasta_mod_file, mod_dict, bounds_dict, ecami_dict, diamond_dict


def cli_prune_seqs():
    parser = argparse.ArgumentParser()
    # fasta_filepath, bounds_file, family, output_folder, source
    parser.add_argument("--fasta", "-a", type=str, help="Fasta file to prune sequences of.")
    parser.add_argument("--bounds", "-b", type=str, help="Bounds file which contains start and stop indices of the "
                                                         "fasta file to prune with.")
    parser.add_argument("--family", "-f", type=str, help="The family to extract sequences of for pruning. The pruned "
                                                         "file will only contain modules from this family.")
    parser.add_argument("--out_folder", "-o", type=str, help="The output folder to write the pruned fasta file and "
                                                             "modified bounds file to. If not given, the fasta data "
                                                             "will be written to standard console output.",
                        default=None)
    parser.add_argument("--source", "-s", type=str, help="The program which generated the bounds file. Currently, only"
                                                         "dbcan hmmer.out and pfam files are supported.",
                        choices=["dbcan", "pfam"])
    args = parser.parse_args()

    pruned, bounds_dict = filter_prune(args.fasta, args.bounds, args.family, args.out_folder, args.source, prune=True,
                                       accession_dict=None)
    if args.out_folder is None:
        out_handle = sys.stdout
    else:
        out_filename = re.sub(r"\.fasta", "_pruned.fasta", os.path.basename(args.fasta))
        out_handle = os.path.join(args.out_folder, out_filename)

    try:
        if len(pruned) > 0:
            SeqIO.write(pruned, out_handle, "fasta")
        else:
            print("# WARNING: No sequences matched the specified family! File not written.")
    except IOError as error:
        print("ERROR: ", error.args[0])
        print(f"ERROR: Cannot write pruned sequences to file {out_handle}.")


if __name__ == "__main__":
    cli_prune_seqs()
    # input_fasta = "PL9_CHARACTERIZED_cazy.fasta"
    # family_test = "PL9"
    # mode_test = Mode.CHARACTERIZED
    # in_path = os.path.join(os.getcwd(), "../../output", family_test, mode_test.name, "cazy", "PL9_CHARACTERIZED_cazy.fasta")
    # out_folder = os.path.join(os.getcwd(), "../../output", family_test, mode_test.name, "dbcan2")
    # # main(in_path, family_test, out_folder, mode_test, )
