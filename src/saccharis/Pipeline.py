###############################################################################
# Main pipeline for SACCHARIS 2.0
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
"""
The Pipeline module contains the code for running the main pipeline. If you want to run the pipeline to create and
optionally render a phylogeny on a single CAZy family with a set of user sequences and settings, you can call the
`single_pipeline` method. The Pipeline module does not include code which processes data, but handles data flow to and
from the code in other modules which handle their respective steps in the pipeline.

There are seven steps in the core pipeline:
    - Step One - CAZyme download
    - Step Two - Load user sequence data
    - Step Three - Module extraction & pruning
    - Step Four - Multiple sequence alignment
    - Step Five - Mutation modelling
    - Step Six - Build Tree
    - Step Seven - Render Tree (optional, requires separately installed R package)

## Step One - CAZyme download
Download specified CAZyme metadata from CAZy, then download sequence data from NCBI. Handled by the `Cazy_Scrape`
module.
## Step Two - Load user sequence data
## Step Three - Module extraction & pruning
## Step Four - Multiple sequence alignment
## Step Five - Mutation modelling
## Step Six - Build Tree
## Step Seven - Render Tree (optional, requires separately installed R package)
"""
import datetime
import json
import logging
import math
import os
import shutil
import sys
import time

from PyQt5.QtCore import pyqtSignal

from saccharis.ParseUserSequences import merge_data_sources
from saccharis.utils.FamilyCategories import Matcher
from saccharis.utils.PipelineErrors import PipelineException, UserError, make_logger
from saccharis import CazyScrape
from saccharis import ChooseAAModel
from saccharis.ExtractAndPruneCAZymes import main as extract_pruned
from saccharis import FastTree_Build
from saccharis import Muscle_Alignment
from saccharis import RAxML_Build
from saccharis.Rendering import render_phylogeny
from saccharis.utils.AdvancedConfig import get_user_settings, get_log_folder, get_version, save_dict_to_file
from saccharis.utils.FamilyCategories import check_deleted_families
from saccharis.utils.Formatting import make_metadata_dict, format_time, CazymeMetadataRecord
from saccharis.CazyScrape import Domain


def single_pipeline(family: str, output_folder: str | os.PathLike,
                    scrape_mode: CazyScrape.Mode = CazyScrape.Mode.CHARACTERIZED, domain_mode: Domain | int = 0b11111,
                    threads: int = math.ceil(os.cpu_count() * 0.75),
                    tree_program: ChooseAAModel.TreeBuilder = ChooseAAModel.TreeBuilder.FASTTREE,
                    get_fragments: bool = False, prune_seqs: bool = True, verbose: bool = False,
                    force_update: bool = False, user_files: list[str | os.PathLike] = None,
                    ncbi_genomes: list[str] = None, ncbi_genes: list[str] = None, auto_rename: bool = False,
                    settings: dict = None, gui_step_signal: pyqtSignal = None,
                    logger: logging.Logger = logging.getLogger(), skip_user_ask: bool =False, render_trees: bool = False,
                    ask_func=None):
    """
    Runs the SACCHARIS pipeline on a single CAZyme family with optional user sequences to create a phylogenetic tree of
    sequences from CAZy.org and user FASTA files.

    :param family: The family whose sequences will be downloaded from http://www.CAZy.org.
    :param output_folder: The folder which final and intermediate results will be saved to.
    :param scrape_mode: A filter to choose what characterization level of cazymes will be downloaded from
    http://www.CAZy.org. Default is characterized. When calling as a function, use the Cazy_Scrape.Mode enum type to
    clearly indicate the scrape mode.
        Allowable modes:
            Cazy_Scrape.Mode.CHARACTERIZED
            Cazy_Scrape.Mode.ALL_CAZYMES
            Cazy_Scrape.Mode.STRUCTURE
    :param domain_mode: A filter to choose which organism domains will have sequences downloaded from the CAZy
    database. Default mode is all domains. When calling as a function, use the `saccharis.Cazy_Scrape.Domain` enum type to specify
    a single domain, or combine them together with the | operator to form desired combinations of domains. This
    argument is functionally just a binary number using a bitmask to indicate whether ot include different domains at
    each position.
        Allowable modes:
        -    Cazy_Scrape.Domain.ARCHAEA
        -    Cazy_Scrape.Domain.BACTERIA
        -    Cazy_Scrape.Domain.EUKARYOTA
        -    Cazy_Scrape.Domain.VIRUSES
        -    Cazy_Scrape.Domain.UNCLASSIFIED

        Any combination of domains can be joined together with | operator like the following examples:
        -    Cazy_Scrape.Domain.ARCHAEA | Cazy_Scrape.Domain.BACTERIA
        -    Cazy_Scrape.Domain.EUKARYOTA | Cazy_Scrape.Domain.VIRUSES | Cazy_Scrape.Domain.BACTERIA


    :param threads: Some tools(e.g. RAxML) allow the use of multi-core processing.  Set a number in here from 1 to
    `<max_cores>`. The default is set at 3/4 of the number of logical cores reported by your operating system. This is
    to prevent lockups if other programs are running.
    :param tree_program: Choice of tree building program. FastTree is the default because it is substantially faster
    than RAxML. RAxML may take days or even weeks to build large trees, but can build higher quality
    trees than FastTree. Both RAxML and RAxML-NG are supported. When calling as a function, use the
    `ChooseAAModel.TreeBuilder` enum type to clearly indicate the tree building program.
        Allowable modes:
            ChooseAAModel.TreeBuilder.RAXML
            ChooseAAModel.TreeBuilder.FASTTREE
            ChooseAAModel.TreeBuilder.RAXML_NG
    :param get_fragments:
    :param prune_seqs:
    :param verbose:
    :param force_update:
    :param user_file:
    :param genbank_genomes:
    :param genbank_genes:
    :param auto_rename:
    :param settings:
    :param gui_step_signal:
    :param merged_dict:
    :param logger:
    :param skip_user_ask:
    :param render_trees:
    :return:
    """
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
    msg = f"Pipeline start, version: {get_version()} "
    logger.debug(msg)

    if verbose:
        logger.setLevel(logging.INFO)

    family = family.upper()
    matcher = Matcher()
    if not matcher.valid_cazy_family(family):
        raise UserError(f"Invalid CAZyme family: {family}")
    check_deleted_families(family)

    # Set up family and group folders
    domain_dir_name = "ALL_DOMAINS" if domain_mode == 0b11111 else \
        ''.join([dom.name[0] for dom in CazyScrape.Domain if dom.value & domain_mode])
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
    #     ''.join([dom.name[0] for dom in CazyScrape.Domain if dom.value & domain_mode])
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
    cazy_file, cazy_metadata, cazy_stats, cazy_seqs = \
        CazyScrape.main(family, cazy_folder, scrape_mode, get_fragments, verbose, force_update, ncbi_query_size,
                         domain_mode, skip_ask=bool(gui_step_signal) or skip_user_ask, logger=logger)
    cazy_t = time.time()
    print("Completed Cazy Extract")
    print("==============================================================================\n")
    print("*********************************************")
    print("* Cazy/NCBI querying takes:")
    print(format_time(cazy_t - start_t))
    print("*********************************************")

    #######################################
    # Step Two - Combine User, Cazy, and NCBI
    #######################################
    user_t = None  # This line to suppress warning from undeclared variable being accessed

    user_folder = os.path.join(domain_folder, "user")
    if not os.path.isdir(user_folder):
        os.mkdir(user_folder, 0o755)

    all_seqs_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}"

    if user_files is not None or ncbi_genomes is not None or ncbi_genes is not None:

        all_seqs, all_metadata, all_seqs_file_path, user_run_id = \
            merge_data_sources(cazy_seqs, cazy_metadata, user_files, ncbi_genomes, ncbi_genes, user_folder,
                               all_seqs_filename, verbose, force_update, auto_rename, logger, skip_user_ask, ask_func)
        user_t = time.time()
        print("Added user FASTA/NCBI sequences to CAZy family sequences")
        print("==============================================================================\n")
        print("*********************************************")
        print("* Appending user sequences takes:")
        print(format_time(user_t - cazy_t))
        print("*********************************************")
    else:  # Source is CAZy only
        all_seqs = cazy_seqs
        all_metadata = cazy_metadata
        all_seqs_file_path = cazy_file
        user_run_id = None

    user_run_insert = f"_UserRun{user_run_id:05d}" if user_run_id is not None else ""

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
    print(f"dbCAN processing of {os.path.split(all_seqs_file_path)[1]} is underway...")
    pruned_list, pruned_file, id_convert_dict, bound_dict, ecami_dict, diamond_dict = \
        extract_pruned(all_seqs_file_path, family, dbcan_folder, scrape_mode, force_update, prune_seqs,
                       threads=threads, hmm_cov=hmm_cov, hmm_eval=hmm_eval)

    try:
        pruned_module_list = [id_convert_dict[seq_record.id] for seq_record in pruned_list]
    except KeyError as err:
        msg = f"pruned_list: {pruned_list}"
        logger.error(msg)
        msg = f"id_convert_dict: {pruned_list}"
        logger.error(msg)
        logger.error(err.args[0])
        raise PipelineException("Error with pruned_list conversion. This may be caused by loading incorrectly "
                                "formatted data, try running the pipeline with --fresh to skip loading partial "
                                "run data.") from err

    metadata_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}{user_run_insert}.json"

    final_metadata_filepath = os.path.join(domain_folder, metadata_filename)

    if not force_update and os.path.isfile(final_metadata_filepath):
        try:
            with open(final_metadata_filepath, 'r', encoding="utf-8") as meta_json:
                cazyme_dict = json.loads(meta_json.read())
                final_metadata_dict = {uid: CazymeMetadataRecord(**record) for uid, record in cazyme_dict.items()}
        except IOError as err:
            logger.debug(err)
            # todo: consider renaming the old JSON file to retain that data and writing the new data when this occurs.
            msg = f"Data from previous run unable to be read! " \
                  f"CazymeMetadataRecordFile with read error: {final_metadata_filepath}\n" \
                  f"Falling back to recalculating fresh CazymeMetadataRecords, but NOT overwriting old ones. If you " \
                  f"want to overwrite the old records, run the pipeline again with the --fresh option."
            logger.error(msg)
            final_metadata_dict = make_metadata_dict(all_metadata, pruned_module_list, bound_dict,
                                                     ecami_dict, diamond_dict, logger=logger)
    else:
        final_metadata_dict = make_metadata_dict(all_metadata, pruned_module_list, bound_dict,
                                                 ecami_dict, diamond_dict, logger=logger)
        try:
            with open(final_metadata_filepath, 'w', encoding="utf-8") as meta_json:
                json.dump(final_metadata_dict, meta_json, default=vars, ensure_ascii=False, indent=4)
        except IOError as error:
            raise UserWarning("Problem writing final module metadata information to file. Make sure you have access "
                              "permissions for your output folder, as this is a common source of write errors of this "
                              "type."
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
    aligned_ren_path, aligned_path, aligned_fasttree = Muscle_Alignment.main(pruned_file, cazyme_module_count,
                                                                             family, scrape_mode,
                                                                             muscle_folder, id_convert_dict,
                                                                             force_update=force_update,
                                                                             user_run_id=user_run_id,
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
                                                scrape_mode, threads, tree_program, force_update, user_run_id,
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
    elif tree_program == ChooseAAModel.TreeBuilder.RAXML_NG:
        print(f"RaxML-NG - Tree building of {os.path.split(aligned_ren_path)[1]} is underway")
        raxml_ng_folder = os.path.join(domain_folder, "raxml_ng")
        if not os.path.isdir(raxml_ng_folder):
            os.mkdir(raxml_ng_folder, 0o755)
        tree_path = RAxML_Build.build_tree_raxml_ng(aligned_ren_path, aa_model, raxml_ng_folder, cazyme_module_count,
                                                    threads, force_update, user_run_id, logger)
    elif tree_program == ChooseAAModel.TreeBuilder.RAXML:
        print(f"RaxML - Tree building of {os.path.split(aligned_ren_path)[1]} is underway")
        raxml_folder = os.path.join(domain_folder, "raxml")
        if not os.path.isdir(raxml_folder):
            os.mkdir(raxml_folder, 0o755)
        tree_path = RAxML_Build.main(aligned_ren_path, aa_model, raxml_folder, raxml_cmd, cazyme_module_count, threads,
                                     force_update, user_run_id, logger)
    else:
        raise PipelineException("Undefined tree construction software specified. This is a bug, it should never "
                                "happen, please report to developer through github or otherwise!")
    print("Completed Building of Tree")
    print("==============================================================================\n")
    tree_t = time.time()
    print("*********************************************")
    print("* Tree building takes:")
    print(format_time(tree_t - prottest_t))
    print("*********************************************")

    # Copy the Best tree and settings file to the group directory
    try:
        tree_fname = f"{family}_{scrape_mode.name}_{domain_dir_name}{user_run_insert}_{tree_program.name}.tree"

        final_tree_path = os.path.join(domain_folder, tree_fname)
        shutil.copyfile(tree_path, final_tree_path)

        settings_filename = f"{family}_{scrape_mode.name}_{domain_dir_name}{user_run_insert}_settings-" \
                            f"{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.json"
        run_settings_path = os.path.join(domain_folder, settings_filename)
        save_dict_to_file(settings, settings_path=run_settings_path)
    except IOError as error:
        raise UserWarning("Problem writing final tree output information to file. Make sure you have access "
                          "permissions for your output folder, as this is a common source of write errors of this type."
                          ) from error

    #######################################
    # Step Seven - Render Tree
    #######################################
    if gui_step_signal:
        # noinspection PyUnresolvedReferences
        gui_step_signal.emit(7)
        if sys.gettrace():
            time.sleep(2)  # this is only active while debugging, for gui testing on already run families

    root = "OUT0000000" if "OUT0000000" in final_metadata_dict else None
    if render_trees:
        print(f"rsaccharis - Tree rendering of {family} is underway")

        render_phylogeny(json_file=final_metadata_filepath, tree_file=final_tree_path, output_folder=domain_folder,
                         root=root)

        print("Completed Rendering of Graphics")
        print("==============================================================================\n")
        render_t = time.time()
        print("*********************************************")
        print("* Rendering took:")
        print(format_time(render_t - tree_t))
        print("*********************************************")

    # Final Benchmark tests
    print("==============================================================================\n")
    print("*********************************************")
    print("* Cazy Pipeline Took in Total:")
    print(format_time(tree_t - start_t))
    print("*********************************************\n")

    print(f"Finished Cazy pipeline analysis for group: {scrape_mode.name} of family {family}")
    print("==============================================================================")
    print("==============================================================================\n")
