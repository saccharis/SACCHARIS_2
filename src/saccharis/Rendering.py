import subprocess
from logging import getLogger, Logger


def render_phylogeny(json_file: str, tree_file: str, output_folder: str, logger: Logger = getLogger()):
    args = ['Rscript',  f'"library(rsaccharis);C_load_and_plot_all("{json_file}", "{tree_file}", "{output_folder}")"']
    try:
        subprocess.run(args, check=True)
        logger.info(f"Successfully rendered phylogenetic trees to folder: {output_folder} ")
    except (subprocess.SubprocessError, subprocess.CalledProcessError) as error:
        logger.debug(error)
        logger.warning("Error running Rscript phylogeny rendering code. Check that rsaccharis is installed in R and "
                       "'Rscript' executable is available on PATH. One some systems 'Rscript' needs to be available on "
                       "the system path, not just user path.")
    except Exception as error:
        logger.error(error)
        logger.error(f"Failed to render phylogenetic trees to output folder: {output_folder}")
        print("Check that rsaccharis is installed in R and 'Rscript' executable is available on PATH.")
