###############################################################################
# Main pipeline for SACCHARIS 2.0
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
import datetime
import json
import logging
import math
import os
import shutil
import sys
import time

from PyQt5.QtCore import pyqtSignal

from saccharis.utils.FamilyCategories import Matcher
from saccharis.utils.PipelineErrors import PipelineException, UserError, make_logger
from saccharis import Cazy_Scrape
from saccharis import ChooseAAModel
from saccharis.ExtractAndPruneCAZymes import main as extract_pruned
from saccharis import FastTree_Build
from saccharis import Muscle_Alignment
from saccharis import Parse_User_Sequences
from saccharis import RAxML_Build
from saccharis.Rendering import render_phylogeny
from saccharis.utils.AdvancedConfig import get_user_settings, get_log_folder
from saccharis.utils.AdvancedConfig import save_to_file
from saccharis.utils.FamilyCategories import check_deleted_families
from saccharis.utils.Formatting import make_metadata_dict, format_time


def single_pipeline(family: str, output_folder: str, scrape_mode: Cazy_Scrape.Mode = Cazy_Scrape.Mode.ALL_CAZYMES,
                    domain_mode: int = 0b11111, threads: int = math.ceil(os.cpu_count() * 0.75),
                    tree_program: ChooseAAModel.TreeBuilder = ChooseAAModel.TreeBuilder.FASTTREE,
                    get_fragments: bool = False, prune_seqs: bool = True, verbose: bool = False,
                    force_update: bool = False, user_file=None, genbank_genomes=None, genbank_genes=None,
                    auto_rename: bool = False, settings: dict = None, gui_step_signal: pyqtSignal = None,
                    merged_dict: dict = None, logger: logging.Logger = None, skip_user_ask=False):

    # todo: remove windows block once WSL support is fully implemented
    if sys.gettrace():
        print("Debug is active")
    # else:
    #     if sys.platform.startswith("win"):
    #         raise UserWarning("Windows support is not yet fully implemented. It should be implemented soon, "
    #                           "for now you can install SACCHARIS through WSL and use the CLI, or possibly run the GUI"
    #                           " through WSL on Windows 11, or use a linux OS to run the GUI.")

    if logger is None:
        logger = make_logger("PipelineLogger", get_log_folder(), "pipeline_logs.txt")

    family = family.upper()
    matcher = Matcher()
    if not matcher.valid_cazy_family(family):
        raise UserError(f"Invalid CAZyme family: {family}")
    check_deleted_families(family)

    # Set up family and group folders
    domain_dir_name = "ALL_DOMAINS" if domain_mode == 0b11111 else \
        ''.join([dom.name[0] for dom in Cazy_Scrape.Domain if dom.value & domain_mode])
    if not prune_seqs:
        domain_dir_name += "_noprune"
    if get_fragments:
        domain_dir_name += "_withfrags"
    output_folder = os.path.abspath(output_folder)  # compatibility line to parse paths correctly on windows
    domain_folder = os.path.join(output_folder, f"{family}_{scrape_mode.name}_{domain_dir_name}")
    cazy_folder = os.path.join(domain_folder, "cazy")
    dbcan_folder = os.path.join(domain_folder, "dbcan2")
    muscle_folder = os.path.join(domain_folder, "muscle")
    prottest_folder = os.path.join(domain_folder, "modeltest")
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder, 0o755)
    if not os.path.isdir(domain_folder):
        os.mkdir(domain_folder, 0o755)
    if not os.path.isdir(cazy_folder):
        os.mkdir(cazy_folder, 0o755)
    if not os.path.isdir(dbcan_folder):
        os.mkdir(dbcan_folder, 0o755)
    if not os.path.isdir(muscle_folder):
        os.mkdir(muscle_folder, 0o755)
    if not os.path.isdir(prottest_folder):
        os.mkdir(prottest_folder, 0o755)
    # # Older code, this block makes a nested folder structure, which was not the preferred structure.
    # # Keeping it around in case we add this as an alternate option.
    # family_folder = os.path.join(output_folder, family)
    # group_folder = os.path.join(output_folder, family, scrape_mode.name)
    # domain_dir_name = "ALL_DOMAINS" if domain_mode == 0b11111 else \
    #     ''.join([dom.name[0] for dom in Cazy_Scrape.Domain if dom.value & domain_mode])
    # if not prune_seqs:
    #     domain_dir_name += "_noprune"
    # if get_fragments:
    #     domain_dir_name += "_withfrags"
    # domain_folder = os.path.join(group_folder, domain_dir_name)
    # cazy_folder = os.path.join(domain_folder, "cazy")
    # dbcan_folder = os.path.join(domain_folder, "dbcan2")
    # muscle_folder = os.path.join(domain_folder, "muscle")
    # prottest_folder = os.path.join(domain_folder, "prottest")
    # if not os.path.isdir(output_folder):
    #     os.mkdir(output_folder, 0o755)
    # if not os.path.isdir(family_folder):
    #     os.mkdir(family_folder, 0o755)
    # if not os.path.isdir(group_folder):
    #     os.mkdir(group_folder, 0o755)
    # if not os.path.isdir(domain_folder):
    #     os.mkdir(domain_folder, 0o755)
    # if not os.path.isdir(cazy_folder):
    #     os.mkdir(cazy_folder, 0o755)
    # if not os.path.isdir(dbcan_folder):
    #     os.mkdir(dbcan_folder, 0o755)
    # if not os.path.isdir(muscle_folder):
    #     os.mkdir(muscle_folder, 0o755)
    # if not os.path.isdir(prottest_folder):
    #     os.mkdir(prottest_folder, 0o755)

    if settings is None:
        settings = get_user_settings()
    ncbi_query_size = settings["genbank_query_size"]
    raxml_cmd = settings["raxml_command"]

    start_t = time.time()
    print("==============================================================================")
    print("==============================================================================")
    print(f"Begin pipeline analysis for group: {scrape_mode.name} of family {family}")

    #######################################
    # Step One - Cazy Extract
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(1)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families
    print(f"Cazy Extract is proceeding for {scrape_mode.name} of family {family}...")
    fasta_file, cazymes, cazy_stats = Cazy_Scrape.main(family, cazy_folder, scrape_mode, get_fragments, verbose,
                                                       force_update, ncbi_query_size, domain_mode,
                                                       skip_ask=bool(gui_step_signal) or skip_user_ask, logger=logger)
    cazy_t = time.time()
    print("Completed Cazy Extract")
    print("==============================================================================\n")
    print("*********************************************")
    print("* Cazy/NCBI querying takes:")
    print(format_time(cazy_t - start_t))
    print("*********************************************")

    #######################################
    # Step Two - Combine User and Cazy
    #######################################
    fasta_with_user_file = ""       #
    user_run_id = None              # These lines to suppress warnings from undeclared variables being accessed
    user_t = None                   #
    if user_file is not None:
        user_folder = os.path.join(domain_folder, "user")
        if not os.path.isdir(user_folder):
            os.mkdir(user_folder, 0o755)
        try:
            fasta_with_user_file, user_count, user_run_id = Parse_User_Sequences.run(user_file, fasta_file, user_folder,
                                                                                     verbose, force_update,
                                                                                     auto_rename or skip_user_ask)
        #     todo: replace this with functions that return seq and cazymemetadatarecord lists to more easily concat
        #       mixtures of family(ies?), genbank genomes/genes, and user seqs
        except UserWarning as error:
            logger.warning(error.args[0])
            answer = input("Would you like to continue anyway, without user sequences in the analysis?")
            if not skip_user_ask and (answer.lower() == "y" or answer.lower() == "yes"):
                print("Continuing...")
                logger.info("Continuing...")
            else:
                print("Exiting...")
                logger.info("Exiting...")
                sys.exit()
        user_t = time.time()
        print("Added user sequences to CAZy family sequences")
        print("==============================================================================\n")
        print("*********************************************")
        print("* Appending user sequences takes:")
        print(format_time(user_t - cazy_t))
        print("*********************************************")
    #######################################
    # Step Three - dbCAN, extract & prune
    #######################################
    # todo: Somewhere during this step, prune seqs below threshold length? Is this before or after pruning cazymes?
    #  Do we even care about this feature?
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(3)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families
    hmm_cov, hmm_eval = settings["hmm_cov"], settings["hmm_eval"]
    if user_file is not None:
        print(f"dbCAN processing of {os.path.split(fasta_with_user_file)[1]} is underway...")
        pruned_list, pruned_file, id_convert_dict, bound_dict, ecami_dict, diamond_dict = \
            extract_pruned(fasta_with_user_file, family, dbcan_folder, scrape_mode, force_update, prune_seqs,
                           threads=threads, hmm_cov=hmm_cov, hmm_eval=hmm_eval)
        metadata_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}_UserRun{user_run_id:05d}.json"
    else:
        print(f"dbCAN processing of {os.path.split(fasta_file)[1]} is underway...")
        pruned_list, pruned_file, id_convert_dict, bound_dict, ecami_dict, diamond_dict = \
            extract_pruned(fasta_file, family, dbcan_folder, scrape_mode, force_update, prune_seqs,
                           threads=threads, hmm_cov=hmm_cov, hmm_eval=hmm_eval)
        metadata_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}.json"

    final_metadata_filepath = os.path.join(domain_folder, metadata_filename)
    final_metadata_dict = make_metadata_dict(cazymes, list(id_convert_dict.values()), bound_dict, merged_dict,
                                             ecami_dict, diamond_dict)
    try:
        with open(final_metadata_filepath, 'w', encoding="utf-8") as meta_json:
            json.dump(final_metadata_dict, meta_json, default=vars, ensure_ascii=False, indent=4)
    except IOError as error:
        raise UserWarning("Problem writing final module metadata information to file. Make sure you have access "
                          "permissions for your output folder, as this is a common source of write errors of this type."
                          ) from error

    cazyme_module_count = len(pruned_list)
    extract_t = time.time()
    print("Completed dbCAN Processing")
    print("==============================================================================\n")
    print(f"\tPruned Sequence Count: {cazyme_module_count}")
    print("*********************************************")
    print("* dbCAN processing takes:")
    print(format_time(extract_t - user_t) if user_t else format_time(extract_t - cazy_t))
    print("*********************************************")

    if cazyme_module_count < 1:
        raise PipelineException("ERROR: Zero valid sequences in pruned output! No point in alignment or phylogeny.")
    if cazyme_module_count < 2:
        raise PipelineException("ERROR: Only one valid sequence in pruned output! No point in alignment or phylogeny.")

    #######################################
    # Step Four - Muscle
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(4)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families
    print(f"Muscle alignment of {os.path.split(pruned_file)[1]} is underway...")
    if user_file:
        aligned_ren_path, aligned_path, aligned_fasttree = Muscle_Alignment.main(pruned_file, cazyme_module_count,
                                                                                 family, scrape_mode,
                                                                                 muscle_folder, id_convert_dict,
                                                                                 force_update=force_update,
                                                                                 user_run_id=user_run_id,
                                                                                 threads=threads)
    else:
        aligned_ren_path, aligned_path, aligned_fasttree = Muscle_Alignment.main(pruned_file, cazyme_module_count,
                                                                                 family, scrape_mode,
                                                                                 muscle_folder, id_convert_dict,
                                                                                 force_update=force_update,
                                                                                 threads=threads)

    muscle_t = time.time()
    print("Completed Muscle Alignment")
    print("==============================================================================\n")
    print("*********************************************")
    print("* Muscle alignment takes:")
    print(format_time(muscle_t - extract_t))
    print("*********************************************")

    #######################################
    # Step Five - Mutation modelling
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(5)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families
    print(f"ModelTest-NG tree modeling of {os.path.split(aligned_path)[1]} is underway\n")
    aa_model = ChooseAAModel.compute_best_model(aligned_path, pruned_list, family, prottest_folder, cazyme_module_count,
                                                scrape_mode, "MF", threads, tree_program, force_update, user_run_id,
                                                use_modelTest=True, logger=logger)
    prottest_t = time.time()
    print("Best model found via ModelTest")
    print("==============================================================================\n")
    print("*********************************************")
    print("* ModelTest mutation modeling takes:")
    print(format_time(prottest_t - muscle_t))
    print("*********************************************")

    #######################################
    # Step Six - Build Tree
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(6)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families
    if tree_program == ChooseAAModel.TreeBuilder.FASTTREE:
        print(f"FastTree - Tree building of {os.path.split(aligned_fasttree)[1]} is underway")
        fasttree_folder = os.path.join(domain_folder, "fasttree")
        if not os.path.isdir(fasttree_folder):
            os.mkdir(fasttree_folder, 0o755)
        tree_path = FastTree_Build.main(aligned_fasttree, aa_model, fasttree_folder, force_update, user_run_id, logger)
    else:  # RAxML
        print(f"RaxML - Tree building of {os.path.split(aligned_ren_path)[1]} is underway")
        raxml_folder = os.path.join(domain_folder, "raxml")
        if not os.path.isdir(raxml_folder):
            os.mkdir(raxml_folder, 0o755)
        tree_path = RAxML_Build.main(aligned_ren_path, aa_model, raxml_folder, raxml_cmd, cazyme_module_count, threads,
                                     user_run_id, logger)
    print("Completed Building of Tree")
    print("==============================================================================\n")
    tree_t = time.time()
    print("*********************************************")
    print("* Tree building takes:")
    print(format_time(tree_t - prottest_t))
    print("*********************************************")

    # Copy the Best tree and settings file to the group directory
    try:
        if user_file is not None:
            tree_fname = f"{family}_{scrape_mode.name}_{domain_dir_name}_UserRun{user_run_id:05d}_{tree_program.name}.tree"
        else:
            tree_fname = f"{family}_{scrape_mode.name}_{domain_dir_name}_{tree_program.name}.tree"
        final_tree_path = os.path.join(domain_folder, tree_fname)
        shutil.copyfile(tree_path, final_tree_path)
        if user_file is not None:
            settings_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}_UserRun{user_run_id:05d}_settings-" \
                                f"{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.json"
        else:
            settings_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}_settings-" \
                                f"{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.json"
        run_settings_path = os.path.join(domain_folder, settings_filename)
        save_to_file(settings, settings_path=run_settings_path)
    except IOError as error:
        raise UserWarning("Problem writing final tree outputs information to file. Make sure you have access "
                          "permissions for your output folder, as this is a common source of write errors of this type."
                          ) from error

    #######################################
    # Step Seven - Build Tree
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(7)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families

    render_phylogeny(json_file=final_metadata_filepath, tree_file=final_tree_path, output_folder=domain_folder)


    # Final Benchmark tests
    print("*********************************************")
    print("* Cazy Pipeline Took in Total:")
    print(format_time(tree_t - start_t))
    print("*********************************************\n")

    print(f"Finished Cazy pipeline analysis for group: {scrape_mode.name} of family {family}")
    print("==============================================================================")
    print("==============================================================================\n")
