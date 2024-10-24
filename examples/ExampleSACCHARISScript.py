from pathlib import Path

from saccharis.Pipeline import single_pipeline
from saccharis.CazyScrape import Mode
from saccharis.CazyScrape import Domain

from inspect import getsourcefile
example_folder = Path(getsourcefile(lambda: 0)).parent
example_output_folder = example_folder / "output"
example_user_file = example_folder / "sample_user_fasta_GH5_GH31_GH95.fasta"

single_pipeline("GH31", example_output_folder, scrape_mode=Mode.CHARACTERIZED,
                domain_mode=Domain.BACTERIA | Domain.ARCHAEA, user_files=[example_user_file],
                render_trees=True, auto_rename=True, skip_user_ask=True)
