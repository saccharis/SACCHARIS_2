"""
# What is SACCHARIS?
Sequence Analysis and Clustering of CarboHydrate Active enzymes for Rapid Informed prediction of
Specificity (SACCHARIS), is a python based pipeline designed to improve functional predictions of uncharacterized
sequences for any CAZyme or CBM family currently maintained on the http://www.CAZy.org website or within user-defined
datasets through phylogenetic methods.

# Who are these docs for?
This documentation describes the API and internal methods of the SACCHARIS codebase for developers who intend to use
SACCHARIS via python scripts or contribute to maintenance and development of the codebase. For more end-user focused
usage that describes how to use the command line interface (CLI) and graphical user interface (GUI), please consult
the wiki on our github: https://github.com/saccharis/SACCHARIS_2/wiki

# How do I contribute to SACCHARIS?
Please use github to submit bugs, feature requests, and reach out to developers and maintainers if you want to
contribute. We would prefer new contributers reach out before sending pull requests, but please make sure any pull
requests pass our automated test suite using github actions in your own repo before submitting them.

# Citation
Fraser, A.S.C. et al. (2024). SACCHARIS v2: Streamlining Prediction of Carbohydrate-Active Enzyme Specificities Within
Large Datasets. In: Lisacek, F. (eds) Protein Bioinformatics. Methods in Molecular Biology, vol 2836. Humana,
New York, NY. https://doi.org/10.1007/978-1-0716-4007-4_16

# API
Please note that the API should NOT be considered stable at this point. In the future, we may reexport key methods and
classes in a more stable API module, but we have not yet done so.

For most programmatic uses, we recommend simply using the methods available in the `saccharis.Pipeline` module and
importing necessary and recommended enum types from the modules dedicated to each step, as referenced in the parameter
descriptions.

# Examples
The following examples show how to call SACCHARIS in a script to run the pipeline in similar ways as users would call
it from the command line or GUI.

## Running SACCHARIS on a single family
TBD

## Running SACCHARIS on multiple families
TBD

## Running SACCHARIS on multiple families defined in a file
TBD

## Running SACCHARIS in exploratory mode
TBD

"""