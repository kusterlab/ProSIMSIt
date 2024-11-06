# ProSIMSIt

Combining SIMSI-Transfer and Oktoberfest rescoring for improved multi-batch TMT results.

Hamood, F., Gabriel, W., Pfeiffer, P., Kuster, B., Wilhelm, M., The, M.; Publication pending.

## Test dataset

For testing ProSIMSIt, we recommend downloading the phosphoproteome raw files from this publication:
Zecha, J., Bayer, F. P. et al., _[Decrypting drug actions and protein modifications by dose- and time-resolved proteomics.](https://www.science.org/doi/10.1126/science.ade3925)_ Science (New York, N.Y.) 2023, 380, 93–101.

PRIDE link for .raw files and MaxQuant search results: https://www.ebi.ac.uk/pride/archive/projects/PXD057211/


## ProSIMSIt Instructions

We recommend installing ProSIMSIt via Poetry. To install Poetry, please follow the instructions on the [Poetry website](https://python-poetry.org/docs/).

### Installation
To install ProSIMSIt, clone the repository and install the package:

```bash
git clone
cd ProSIMSIt
poetry install
```

### Running ProSIMSIt
ProSIMSIt can be executed via the command line in the generated environment:

```bash
poetry python -m prosimsit -c /path/to/config.toml
```
