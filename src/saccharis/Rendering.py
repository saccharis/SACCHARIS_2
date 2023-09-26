import subprocess
import sys
from logging import getLogger, Logger


def render_phylogeny(json_file: str, tree_file: str, output_folder: str, logger: Logger = getLogger(),
                     root: str = None):
    print("1")
    try:
        result = subprocess.run(["Rscript", "--version"], check=True)
        print(result.stdout)
        print("2")
    except (subprocess.SubprocessError, subprocess.CalledProcessError) as error:
        logger.debug(error)
        logger.warning("Rscript version command failed")
    json_file_double_slash = json_file.replace('\\', '\\\\')
    tree_file_double_slash = tree_file.replace('\\', '\\\\')
    output_folder_double_slash = output_folder.replace('\\', '\\\\')
    root_arg = f", \'{root}\'" if root else ''
    load_call = f"C_load_and_plot_all(\'{json_file_double_slash}\', \'{tree_file_double_slash}\', " \
                f"\'{output_folder_double_slash}\'{root_arg})"
    args = ['Rscript',  '-e', f'"library(rsaccharis)"', '-e',  f'"{load_call}"']
    try:
        # the run call doesn't work with args as a list because of weird unquoting behaviour of
        # subprocess.run() on windows
        if sys.platform.startswith("win"):
            subprocess.run(' '.join(args), check=True)
        else:
            subprocess.run(args, check=True)
        logger.info(f"Successfully rendered phylogenetic trees to folder: {output_folder} ")
    except (subprocess.SubprocessError, subprocess.CalledProcessError) as error:
        logger.debug(error)
        logger.warning("Error running Rscript phylogeny rendering code. Check that rsaccharis is installed in R and "
                       "'Rscript' executable is available on PATH. One some systems 'Rscript' needs to be available on "
                       "the system path, not just user path.\n"
                       "This does not affect the creation of the pipeline output files, you can still run the "
                       "rsaccharis rendering scripts manually.")
    except Exception as error:
        logger.error(error)
        logger.error(f"Failed to render phylogenetic trees to output folder: {output_folder}")
        print("Check that rsaccharis is installed in R and 'Rscript' executable is available on PATH.")
