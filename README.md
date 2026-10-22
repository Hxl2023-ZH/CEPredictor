# CEPredictor

**CEPredictor: A Machine Learning-Based Tool for Predicting Microbial Chain Elongation Potential**

**Version:** 1.0.0
**Release date:** October 22, 2026
**Language:** R
**Architecture:** Object-oriented programming (R6) with a command-line interface (CLI)

## Overview

CEPredictor is an R-based bioinformatics tool designed to evaluate the potential for microbial chain elongation (CE) using genome-derived functional annotations.

The tool provides a computational framework for identifying and prioritizing microbial genomes, including metagenome-assembled genomes (MAGs), with putative chain elongation capabilities.

CEPredictor supports two primary workflows:

- **Model training:** Train a prediction model using an annotated genome dataset.
- **Model prediction:** Apply a previously trained model to new genome annotation datasets.

The tool is implemented using the R6 object-oriented programming framework and provides a command-line interface for reproducible analyses.

## Installation

Download or clone the CEPredictor repository and navigate to the project directory.

```
git clone <repository-url>
cd CEPredictor
```

CEPredictor requires an R environment with the dependencies used by `CEPredictor.R` installed.

The exact package dependencies and supported R versions should be checked against the released script.

## Usage

### 1. Model Training

Train a prediction model using a genome annotation dataset.

```shell
Rscript CEPredictor.R --input merged_annotation.txt --output train_results --train
```

**Arguments**

| Argument   | Description                              |
| ---------- | ---------------------------------------- |
| `--input`  | Path to the input genome annotation file |
| `--output` | Directory for training results           |
| `--train`  | Enable model training mode               |

The training workflow uses the supplied annotation dataset to construct a prediction model.

The resulting trained model can subsequently be used for predicting chain elongation potential in new genomes.

### 2. Prediction Using a Pretrained Model

Apply an existing model to a new genome annotation dataset.

```shell
Rscript CEPredictor.R --input MAG_annotation_test.txt --output test_results --model train_results/models
```

**Arguments**

| Argument   | Description                                          |
| ---------- | ---------------------------------------------------- |
| `--input`  | Path to the genome annotation dataset to be analyzed |
| `--output` | Directory for prediction results                     |
| `--model`  | Path to the pretrained model directory               |

**Note:** The model directory must contain the model files required by the prediction workflow. The input annotations should follow the feature definitions and formatting requirements used during training.

## Input Data

CEPredictor accepts genome-derived functional annotation data.

Example input filenames:

```
merged_annotation.txt
MAG_annotation_test.txt
```

The required table structure, feature names, annotation conventions, and missing-value handling should be documented according to the implementation in `CEPredictor.R`.

## Prediction Framework

CEPredictor is intended to support genome-based evaluation of microbial chain elongation potential.

Chain elongation is commonly associated with reverse β-oxidation (rBOX), which enables successive two-carbon extensions of acyl-CoA intermediates.

Functional genomic evidence can provide useful information about the potential presence of this metabolic capability. However, the detection of pathway-associated genes alone does not establish experimentally verified chain elongation activity or determine the chain length of the final products.

Accordingly, CEPredictor predictions should be interpreted as **genome-based estimates of metabolic potential**, rather than direct measurements of chain elongation activity.

## Recommended Workflow

```
Genome sequences / MAGs
          |
          v
Functional gene annotation
          |
          v
Genome annotation matrix
          |
          v
     CEPredictor
          |
          +------------------------+
          |                        |
          v                        v
    Model training           Model prediction
          |                        |
          v                        |
    Trained model -----------------+
                                   |
                                   v
                           Prediction results
```

## Example Project Structure

```
tools/CEPredictor/
	├── CEPredictor.R
	├── README.md
	├── merged_annotation.txt
	├── MAG_annotation_test.txt
	├── train_results/
	│   └── models/
	└── test_results/
```

This structure illustrates the command-line workflow and is not an exhaustive listing of automatically generated files.

## Interpretation and Limitations

CEPredictor is designed to prioritize genomes with putative chain elongation potential for further investigation.

Several considerations are important when interpreting predictions:

1. Genomic potential does not necessarily indicate active metabolism under specific environmental or cultivation conditions.
2. Reverse β-oxidation pathways are also involved in short-chain fatty acid production, including butyrate formation.
3. The quality and completeness of genome assemblies may influence the detection of metabolic genes.
4. Experimental validation, including fermentation assays and fatty acid product profiling, is recommended to confirm medium-chain fatty acid production.

## Version History

### v1.0.0 — 2026-10-22

- Initial version of CEPredictor.
- R6-based object-oriented implementation.
- Command-line interface for model training and prediction.
- Support for prediction using previously trained models.

## Citation

Citation information will be provided upon publication or release of the associated manuscript.

## License

License information to be specified by the project maintainers.