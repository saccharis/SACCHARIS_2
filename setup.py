from setuptools import setup
# import yaml
# with open('conda-recipe/meta.yaml', 'r') as file:
#     meta_yaml_data = yaml.safe_load(file)

# todo: match all relevant data labels to the standard metadata format
#  https://packaging.python.org/en/latest/specifications/core-metadata/

setup(name='saccharis',
      # Broadly speaking, we will follow a major.minor.patch.dev semantic versioning structure, where the major release
      # corresponds to major publications.
      # See https://peps.python.org/pep-0440/ for info on what versioning conventions python packages typically follow.
      # Version num info:
      #     -   First digit is a publication release related to a specific academic paper. This will be Saccharis 2 at
      #         the time of writing.
      #     -   Second digit is major public releases which contain major new features, breaking changes that may affect
      #         scripts calling SACCHARIS, etc.
      #     -   Third digit is minor release number for public bugfix or minor feature releases which do not change the
      #         public interface. Any release incrementing the third digit should be fully backwards compatible with
      #         previous releases with the same first and second digits.
      #     -   Fourth digit comes after "dev" and is an internal revision number for testing of new
      #         releases.
      # EXAMPLES:
      # - The third internal version for testing leading up to the first public release could be "2.0.0.dev3"
      # - first public release might be "2.0.0"
      # - First bugfix to public release might be "2.0.1"
      # - The first public release after 2.0.* that breaks anything reliant on 2.0.* will be 2.1.0

      version="2.0.0.dev19",
      build=8,
      # version=meta_yaml_data["package"]["version"],
      # build=meta_yaml_data["build"]["number"],
      description='Bioinformatics tool for automated CAZyme phylogeny construction',
      long_description="This is SACCHARIS 2, a bioinformatics tool for using phylogenetic inference to infer CAZyme "
                       "functionality in genetic sequences.",
      author='Alexander Fraser',
      author_email='alexander.fraser@alumni.ubc.ca',
      url='https://github.com/saccharis/SACCHARIS_2',
      classifiers=[
          # Uncomment the release status below which best describes the state of this release.
          'Development Status :: 3 - Alpha',
          # 'Development Status :: 4 - Beta',
          # 'Development Status :: 5 - Stable',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
      ],
      include_package_data=True,
      entry_points={
          "console_scripts": [
              "saccharis = saccharis.CLI:cli_main",
              "saccharis.make_family_files = saccharis.utils.FamilyCategories:cli_main",
              "saccharis.add_family_category = saccharis.utils.FamilyCategories:cli_append_user_family",
              "saccharis.rename_user_file = saccharis.utils.UserFastaRename:cli_main",
              "saccharis.prune_seqs = saccharis.ExtractAndPruneCAZymes:cli_prune_seqs",
              "saccharis.screen_cazome = saccharis.ScreenUserFile:cli_cazome",
              "saccharis.show_family_categories = saccharis.utils.FamilyCategories:show_categories",
              "saccharis-gui = saccharis.gui.PipelineGUI:main",
              "saccharis.config = saccharis.utils.AdvancedConfig:cli_config"
            ]
      },
      package_data={
          # If any package or subpackage contains *.txt or *.yaml files, include
          # them:
          "": ["*.png", "*.yaml"],
          # # Include any *.msg files found in the "hello" package (but not in its
          # # subpackages):
          # "hello": ["*.msg"],
          # # Include any *.csv files found in the "hello.utils" package:
          # "hello.utils": ["*.csv"],
          # # Include any *.dat files found in the "data" subdirectory of the
          # # "mypkg" package:
          # "mypkg": ["data/*.dat"],
      },
      license='GPL v3',
      install_requires=[
          'beautifulsoup4',
          'biopython',
          'dbcan',
          'lxml',
          'ncbi-datasets-pylib',
          'python-dotenv',
          'wget',
          'requests',
          'setuptools',
          'psutil',
          'pyqt5',
          'PyQt5-sip',
      ],
      python_requires='>=3.8',
      zip_safe=False
      )
