# BrainRep
This script codes for a predictive algorithm that uses gene expression data of human brain cells and regulatory network information to predict transcription factors for inducing cell reprogramming between specific brain cells and neuronal subtypes.

## Files
- **brain_cells_data.json** : file containing the identifiers of the brain cells and neuronal subtypes that can be used as source cell or desired cell --types, and the brain region where they are present <br>
- **brain_rep.py** : main script <br>
- **br_utils.py** : script containing helper functions <br>
- **brain_qc.h5ad** : raw gene expression file (NOT uploaded to GitHub) <br>
- **trimmed_means.csv** : log-transformed gene expression file (NOT uploaded to GitHub) <br>

## Usage
brain_rep.py [-h] -s SOURCECELL -d DESIREDCELL [-o ORGANISM]

## Example
python3 brain_rep.py -s "Neuron 001" -d "Neuron 051"
