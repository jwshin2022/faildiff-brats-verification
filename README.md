# FailDiff-BraTS numerical verification

Standalone verification and aggregation software for **Error transfer and certification saturation in a brain tumor segmentation correction system** by Jee-Woo Shin and Chung-Pyo Hong. The manuscript is under author review; this repository does not imply acceptance or publication.

## Scope

This repository contains code only. It contains no medical images, annotation arrays, record-level data, private identifiers, original input hashes, model weights, credentials, or original segmentation training implementation. The MIT license covers the software here, not third-party data or model assets.

The verifier reconstructs regional Dice, voxel transitions, retrospective singleton outcomes, interval saturation, an algebraic joint-exclusion condition, 18 arithmetic diagnostic conditions, and exploratory group bootstrap intervals from integer summary inputs. It does not rerun image preprocessing, model training, inference, or numerical optimization traces.

## Installation and public synthetic tests

Tested with Python 3.12 and NumPy 2.3.5. Create an isolated environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify_results.py --synthetic --out synthetic_result.json
python test_revision_analysis.py
```

Use a new output path for each execution. The first test uses an artificial group. The second exhaustively compares direct masks and disjoint-count aggregation in 4,096 configurations, including empty regions. Neither requires study data.

## Reproducing study numbers

The record-level numerical inputs are not in this public repository. Request the supporting numerical derivatives from the corresponding author, Chung-Pyo Hong (cphong@hoseo.edu), under applicable data-use and institutional requirements. Public code alone is not full study-data access. See DATA_SCHEMA.md.

```sh
python verify_results.py --data /approved/path/minimum_data.json --out numerical_results.json
python revision_analysis.py --data /approved/path/minimum_data.json --out new_analysis_directory
```

The revision output directory must not already exist. It contains group-level plotting statistics and provenance, which must remain in authorized storage and must not be uploaded by default. The script reports every threshold (0.5, 0.7, 0.9); it does not select a favorable threshold.

Primary bootstrap: 20,000 group draws, PCG64 seed 20261008. Worst-B exclusions: post hoc sensitivity, seed 20261103. Groups are source/image units, not independently verified patients. Numerical reproduction does not resolve source attribution, prove independent sampling, or establish clinical safety. Counts are exact; floating-point bootstrap values should be compared using documented tolerances and dependency versions.

## Development and license

MIT License. Research model outputs and original source bundles remain unchanged. Software and access questions should be directed to the corresponding author. No patient or medical data should be attached to public issues.
