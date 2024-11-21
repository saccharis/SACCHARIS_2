###############################################################################
# Main command line interface for SACCHARIS 2.0
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
import argparse
import math
import os
import sys

from saccharis.Cazy_Scrape import Mode, Domain
from saccharis.ChooseAAModel import TreeBuilder
from saccharis.ParseUserSequences import concatenate_multiple_fasta
from saccharis.Pipeline import single_pipeline
from saccharis.ScreenUserFile import choose_families_from_fasta
from saccharis.utils.AdvancedConfig import MultilineFormatter, get_log_folder, get_version
from saccharis.utils.FamilyCategories import Matcher, get_category_list, load_family_list_from_file
from saccharis.utils.PipelineErrors import UserError, PipelineException, NewUserFile, make_logger
from saccharis.utils.Formatting import rename_metadata_dict_ids


def cli_main():
    logger = make_logger("CLILogger", get_log_folder(), "cli_logs.txt")
    parser = argparse.ArgumentParser(description=f'SACCHARIS version: {get_version()}\n'
                                                 f' SACCHARIS 2 is a tool to analyze CAZyme families.',
                                     formatter_class=MultilineFormatter,
                                     epilog="The following list is of additional utilities and commands available to "
                                            "interact with saccharis. Most commands (not saccharis-gui) have their own "
                                            "help menu, similarly accessed via \"<command_name> -h\" or "
                                            "\"<command_name> --help\", e.g. \"saccharis.make_family_files -h\"|n "
                                            "COMMAND LIST:|n "
                                            "-\tsaccharis.make_family_files|n "
                                            "-\tsaccharis.add_family_category|n "
                                            "-\tsaccharis.rename_user_file|n "
                                            "-\tsaccharis.prune_seqs|n "
                                            "-\tsaccharis.screen_cazome|n "
                                            "-\tsaccharis.show_family_categories|n "
                                            "-\tsaccharis.config|n "
                                            "-\tsaccharis.update_db|n "
                                            "-\tsaccharis-gui"
                                     )
    parser.add_argument('--version', "-v", action='version', version=f"SACCHARIS {get_version()}")
    parser.add_argument('--directory', "-o", type=str, default=os.path.join(os.getcwd(), "output"), help='You can set '
                        'a predefined output directory with this flag, either a full path or a subfolder of the CWD.  '
                        'Default is <Current Working Directory (CWD)>/output. If you specify an absolute file path the '
                        'end directory will be used. If you specify a relative file path(e.g. just a folder name), it'
                        ' will be a subdirectory of the CWD.')
    family_group = parser.add_mutually_exclusive_group(required=True)
    family_group.add_argument("--family", "-f", type=str, help="This is a single family name.\n"
                              "-> eg. \"GH43\". Cannot use with --family_file, --family_category, or explore")
    family_group.add_argument("--family_file", type=argparse.FileType('r'), help="This is a file containing a list of "
                              "families you would like to run the pipeline on sequentially. Cannot use with "
                              "--family/-f, --family_category, or explore.")
    family_group.add_argument("--family_category", type=str, help="This accepts the name of a list of families "
                              "contained in the \"family_categories.json\" config file. Custom groupings can be added"
                              " to this file by the user via text editor, provided that you follow the same "
                              "formatting. The formatting is standard json format as used by the default python json"
                              " encoder for a dict of lists of strings. Cannot use with "
                              "--family/-f, --family_file, or explore.")
    family_group.add_argument("--explore", "-x", action="store_true", help="This is a boolean flag which enables an "
                              "exploratory mode to screen the user sequence file for cayme families and then ask the "
                              "user which of the families in their sequence file to run the pipeline on.Cannot use "
                              "with --family/-f, --family_category, or --family_file.")
    parser.add_argument("--subfamily", type=int, help="This is a subfamily number, which cazy represents "
                        "as a number appended to the family, e.g. \"--family GH43 --subfamily 1\" accesses"
                        " the family which the CAZy database identifies as GH43_1. Optional argument, "
                        "can only use with --family, not --family_file, --family_category, or --explore.",
                        default=None)
    parser.add_argument("--cazyme_mode", "-c", type=str, help="This is the characterization level of cazymes you wish "
                                                              "to query from the CAZy database. Allowable modes:\n"
                                                              "\t characterized\n"
                                                              "\t structure\n"
                                                              "\t all_cazymes\n",
                        choices=["characterized", "structure", "all_cazymes"], default="characterized")
    parser.add_argument("--domain", "-d", type=str, help="This is the domain of the organisms whose cazymes you wish "
                                                         "to query from the CAZy database. Default mode is all domains."
                                                         "\nAllowable modes:\n"
                                                         "\t archaea\n"
                                                         "\t bacteria\n"
                                                         "\t eukaryota\n"
                                                         "\t viruses\n"
                                                         "\t unclassified\n"
                                                         "\t all\n"
                                                         "You can specify as many of these domain options as you wish "
                                                         "by separating them with a space.\n"
                                                         "\t e.g. '-d archaea bacteria' is a valid domain argument "
                                                         "which will include cazyme sequences from organisms in both "
                                                         "the archaea and bacteria domains.",
                        nargs='+', default="all", choices=["archaea", "bacteria", "eukaryota", "viruses",
                                                           "unclassified", "all"])
    parser.add_argument("--seqfile", "-s", nargs='+', type=str, help="If you would like to add your own sequences to "
                        "this run - this is your chance.  Sequences MUST be in FASTA FORMAT - if they are not the "
                        "script will fail.  Make sure to include path with filename. Multiple FASTA files are "
                        "supported and will be merged together, with metadata of source files saved for future use.")
    # genbank_group = parser.add_argument_group(required=False)
    # genbank_group.add_argument("--ncbi_genome", "-g", nargs='+', type=str, help="If you would like to add sequences "
    #                            "from a genbank genome, this is the argument to do so. Just add the genome assembly ID "
    #                            "(these start with GCA_ or GCF_) and the protein coding sequences will be downloaded in "
    #                            "FASTA format and added as user sequences to your results. Multiple genomes are "
    #                            "supported, each id separated with a space, and will be merged together with all the "
    #                            "user sequences into the resultant screening/tree, with metadata of source saved for "
    #                            "future use.")
    # todo: implement gene fasta download
    # genbank_group.add_argument("--ncbi_genes", "-n", nargs='+', type=str, help="If you would like to add sequences "
    #                            "from a list of genbank gene IDs, this is the argument to do so. Just add the GenBank "
    #                            "gene IDs as a space separated list to this argument and the protein sequences will be "
    #                            "downloaded in FASTA format and merged together with all the user sequences into the "
    #                            "resultant screening/tree, with metadata of source saved for future use.")
    parser.add_argument("--rename_user", "-u", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means the program will not automatically rename user sequence headers "
                        "to conform with the user sequence ID format. When this argument is included, this will occur "
                        "automatically, without prompting. When not included, the program will prompt the user if they "
                        "wish to prepend their FASTA headers with the correct ID format.")
    parser.add_argument("--fresh", "-n", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means existing data will be reused to speed up analysis. When included,"
                        " this options forces data to be redownloaded from CAZy and"
                        " NCBI even if present and all analyses to be performed again with fresh data. Saved data "
                        "from previous runs with the same family, cazyme mode, and domain settings will be deleted and "
                        "overwritten.")
    parser.add_argument("--verbose", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means verbose output is hidden. If you would like verbose output "
                        "(particularly useful to explore when certain sequences from CAZy are not being included"
                        " in an analysis), include this flag in your call.")
    parser.add_argument("--skip_prune", "-k", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means the sequences will be pruned to CAZyme boundaries. If you would"
                        "like to skip the pruning step and use full genbank entries for alignment and tree building, "
                        "include this flag in your call.")
    parser.add_argument("--fragments", "-m", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means fragments are left out by default. If you would like to include "
                        "fragment sequences from CAZY, include this flag in your call.")
    parser.add_argument("--threads", "-t", type=int, default=math.ceil(os.cpu_count()*0.75),
                        help="Some tools(e.g. RAxML) allow the use of multi-core processing.  Set a number in here from"
                             " 1 to <max_cores>. The default is set at 3/4 of the number of logical cores reported by "
                             "your operating system. This is to prevent lockups if other programs are running.",
                        choices=range(1, os.cpu_count()+1))
    parser.add_argument("--tree", "-e", default="fasttree", choices=["fasttree", "raxml", "raxml_ng"], help="Choice of "
                        "tree building program. FastTree is the default because it is substantially faster than RAxML. "
                        "RAxML may take days or even weeks to build large trees, but sometimes builds slightly higher "
                        "quality trees than FastTree. Both RAxML and RAxML-NG are supported.")
    parser.add_argument("--skip_user_ask", action="store_true", help="This is a boolean flag that by default is set to "
                        "False which when true skips querying the user for input. If a question would be posed to the "
                        "user, the question will be skipped and a reasonable default action will occur if possible or "
                        "the program will exit if necessary. An example is that normally when no NCBI API key is "
                        "detected the program prompts the user for their NCBI API key. If this setting is true, the "
                        "program instead continues without an API key, slowing down queries but otherwise functioning. "
                        "It is recommended to use this option when running SACCHARIS on a cluster or in automated "
                        "fashion to prevent failed runs in environments where querying the user is not possible.")
    parser.add_argument("--render", "-r", action="store_true", help="This is a boolean flag that by default is set to "
                        "False which when true automatically tries to render phylogenetic trees using the rsaccharis R "
                        "library. The rsaccharis package must be installed and the rscript command available on the "
                        "PATH variable for this function to work, see https://github.com/saccharis/rsaccharis for "
                        "details.")

    args = parser.parse_args()

    # validate args
    cazyme_mode = Mode[args.cazyme_mode.upper()]
    fragments = args.fragments
    verbose_arg = args.verbose
    refresh = args.fresh
    skip_user_ask = args.skip_user_ask
    num_threads = args.threads if args.threads <= os.cpu_count() else os.cpu_count()
    render_trees = args.render

    domain_val = 0b0
    if isinstance(args.domain, str):
        if args.domain.upper() == "ALL":
            domain_val = 0b11111
        else:
            domain_val |= Domain[args.domain.upper()].value
    else:
        for item in args.domain:
            if item.upper() == "ALL":
                domain_val = 0b11111
                break
            domain_val |= Domain[item.upper()].value
    prune = not args.skip_prune
    tree_prog = TreeBuilder[args.tree.upper()]
    rename = args.rename_user

    # set output path
    if args.family or args.explore:
        base_dir = os.path.abspath(args.directory)
        output_path = os.path.abspath(args.directory)
    elif args.family_category:
        base_dir = os.path.abspath(args.directory)
        output_path = os.path.join(base_dir, args.family_category)
    elif args.family_file:
        base_dir = os.path.abspath(args.directory)
        output_path = os.path.join(base_dir, os.path.splitext(os.path.basename(args.family_file.name))[0])
    else:
        raise Exception("Something has gone wrong with the command line input parsing while reading output path info.")

    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # user_path = os.path.abspath(args.seqfile) if args.seqfile else None
    if args.seqfile is None:
        user_path = None
        user_merged_dict = None
    elif type(args.seqfile) == list:
        # todo: don't handle genbank fasta download here for both genome and genes, delete the code that did this
        #  in the following function call because it's not easily extensible
        user_path, user_merged_dict, user_seqs = concatenate_multiple_fasta(args.seqfile, output_folder=output_path)
    else:
        raise Exception("Error parsing user sequence file(s) from command line. This shouldn't happen, "
                        "please report as a bug through github.")

    # if args.ncbi_genome is not None:
    #     ncbi_genome = args.ncbi_genome
    # #     todo: this probably shouldn't even be here tbh
    #     # Download seqs from NCBI for given genomes
    #     genome_storage_dir = os.path.join(os.path.expanduser('~'), "saccharis", "ncbi_downloads")
    #     genome_seqs, genome_source = download_proteins_from_genomes(ncbi_genome, out_dir=genome_storage_dir, logger=logger)
    #     # origin_dict += genome_source
    #     # all_seqs += genome_seqs

    matcher = Matcher()
    if args.family:
        family_arg = args.family.upper()
        if not matcher.valid_cazy_family(family_arg):
            print(f"ERROR: Invalid family argument: \"{family_arg}\" \n"
                  f"Please input a valid family: PL*, GH*, GT*, CE*, or AA*, where * is a number.")
            sys.exit(3)
        elif family_arg.__contains__("_") and args.subfamily:
            print(f"ERROR: Family argument \"{family_arg}\" seems to already contain a subfamily and you have specified"
                  f" subfamily \"{str(args.subfamily)}\". Please use ONLY ONE of these methods to specify a subfamily.")
            sys.exit(3)
        if args.subfamily:
            family_arg += f"_{args.subfamily}"

    elif args.family_category:
        if args.subfamily:
            print("ERROR: Cannot use subfamily argument with family categories. Instead, you can customize the JSON "
                  "file which contains the family categories and add subfamilies to default or custom lists using "
                  "\"<family>_<subfamily>\" syntax. \n"
                  "\te.g. subfamily 1 of GH43 is \"GH43_1\" and can be added to category lists specifically")
            sys.exit(3)
        try:
            fam_list = get_category_list(args.family_category)
        except UserError as error:
            print(error.msg)
            sys.exit(3)

        for fam in fam_list:
            if not matcher.valid_cazy_family(fam):
                print(f"ERROR: Invalid family argument read from family category file: \"{fam}\"\n"
                      f"\tPlease input a valid family: PL*, GH*, GT*, CE*, or "
                      f"AA*, where * is a number.")
                sys.exit(3)

    elif args.explore:
        if args.subfamily:
            print("ERROR: Cannot use subfamily argument with explore. Instead, you can select the family_subfamily "
                  "categories found in the file to run the pipeline on to make trees.")
            sys.exit(3)
        if not args.seqfile:
            print("ERROR: Cannot run exploratory mode without a user sequence file!")
            print("Exiting...")
            sys.exit(3)
        try:
            fam_list = choose_families_from_fasta(user_path, output_path, num_threads)
        except PipelineException as error:
            print(f"ERROR: {error.msg}")
            print("Exiting...")
            sys.exit(2)

    elif args.family_file:  # family file
        if args.subfamily:
            print("ERROR: Cannot use subfamily argument with family files. Instead, you can edit your input file which "
                  "contains the family categories and add subfamilies using \"<family>_<subfamily>\" syntax. \n"
                  "\te.g. subfamily 1 of GH43 is \"GH43_1\" and can be added to category lists in that manner")
            sys.exit(3)

        try:
            fam_list = load_family_list_from_file(args.family_file)
        except IOError as error:
            print(f"ERROR: {error.msg}")
            print("ERROR: Error loading data from family file.")
            sys.exit(3)
        except UserError as error:
            print(f"ERROR: {error.msg}")
            print("Exiting...")
            sys.exit(3)
    else:
        raise Exception("Something has gone wrong with command line input parsing while reading family information.")

    if args.family:
        # todo: Refactor this section to only have one single_pipeline call. This whole section of try and excepts is
        #  awful, normal flow control using the NewUserFile exception was a bad idea and bad practice. Single family
        #  should just go into the fam_list and user file control flow should not be exception based. But this
        #  refactoring will take some time and is low priority right now. Probably have to make some kind of wrapper
        #  function for the pipeline to get rid of the exception control flow or something else tbd.
        try:
            single_pipeline(family_arg, output_path, cazyme_mode, domain_mode=domain_val, threads=num_threads,
                            tree_program=tree_prog, get_fragments=fragments, prune_seqs=prune, verbose=verbose_arg,
                            force_update=refresh, user_file=user_path, auto_rename=rename, merged_dict=user_merged_dict,
                            logger=logger, skip_user_ask=skip_user_ask, render_trees=render_trees)
        except NewUserFile as file_msg:
            user_path = file_msg.msg
            user_merged_dict = rename_metadata_dict_ids(user_path, user_merged_dict)
            try:
                single_pipeline(family_arg, output_path, cazyme_mode, domain_mode=domain_val, threads=num_threads,
                                tree_program=tree_prog, get_fragments=fragments, prune_seqs=prune, verbose=verbose_arg,
                                force_update=refresh, user_file=user_path, auto_rename=rename,
                                merged_dict=user_merged_dict, logger=logger, skip_user_ask=skip_user_ask,
                                render_trees=render_trees)
            except PipelineException as pipe_error:
                logger.error(pipe_error.msg)
                logger.debug(pipe_error.__traceback__)
                logger.error(f"Something went wrong running the SACCHARIS pipeline on family: {family_arg}")
        except PipelineException as pipe_error:
            logger.error(pipe_error.msg)
            logger.debug(pipe_error.__traceback__)
            logger.error(f"Something went wrong running the SACCHARIS pipeline on family: {family_arg}")
    else:
        print("Beginning multiple pipeline runs for each of the following families:", fam_list)
        for family_arg in fam_list:
            try:
                single_pipeline(family_arg, output_path, cazyme_mode, domain_mode=domain_val, threads=num_threads,
                                tree_program=tree_prog, get_fragments=fragments, prune_seqs=prune, verbose=verbose_arg,
                                force_update=refresh, user_file=user_path, auto_rename=rename,
                                merged_dict=user_merged_dict, logger=logger, skip_user_ask=skip_user_ask,
                                render_trees=render_trees)
            except NewUserFile as file_msg:
                user_path = file_msg.msg
                user_merged_dict = rename_metadata_dict_ids(user_path, user_merged_dict)
                try:
                    single_pipeline(family_arg, output_path, cazyme_mode, domain_mode=domain_val, threads=num_threads,
                                    tree_program=tree_prog, get_fragments=fragments, prune_seqs=prune,
                                    verbose=verbose_arg, force_update=refresh, user_file=user_path, auto_rename=rename,
                                    merged_dict=user_merged_dict, logger=logger, skip_user_ask=skip_user_ask,
                                    render_trees=render_trees)
                except PipelineException as pipe_error:
                    logger.error(pipe_error.msg)
                    logger.debug(pipe_error.__traceback__)
                    logger.error(f"Something went wrong running the SACCHARIS pipeline on family: {family_arg}")
                    print("\t Continuing to run SACCHARIS pipeline on remaining families...")
                    logger.info("\t Continuing to run SACCHARIS pipeline on remaining families...")
            except PipelineException as pipe_error:
                logger.error(pipe_error.msg)
                logger.debug(pipe_error.__traceback__)
                logger.error(f"Something went wrong running the SACCHARIS pipeline on family: {family_arg}")
                print("\t Continuing to run SACCHARIS pipeline on remaining families...")
                logger.info("\t Continuing to run SACCHARIS pipeline on remaining families...")


if __name__ == "__main__":
    cli_main()
