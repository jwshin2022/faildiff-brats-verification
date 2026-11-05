# Additional distribution and positive-score analyses

This extension uses the same restricted integer-statistic schema. It reconstructs composite regional prediction volumes and true-positive overlaps, then computes Dice. Singleton Dice differences are never added. Only disjoint-cell summaries are supported.

```bash
python test_policy_sensitivity_analysis.py
python policy_sensitivity_analysis.py --data /approved/location/minimum_data.json --root /new/local/output
```

The exhaustive test compares 12,544 selected-subset calculations with direct two-voxel masks, including empty regions. Actual-data execution verifies all 93 stored B count/transition results before evaluating confidence thresholds 0.5, 0.7 and 0.9 and the positive-score-only B sensitivity. All four policies are reported.

Bootstrap: 20,000 group resamples, PCG64, seed 20261008 for complete sets; 20261104 for post hoc within-policy tail exclusions. Original B/point_B tail analyses in distribution_analysis.py retain seed 20261103 and their original shared exclusion set.

The output private/ folder contains individual derived measurements and binding hashes. Keep it internal. The analysis/ output contains aggregate summaries and anonymized extreme-outcome summaries; external sharing still requires appropriate review. No medical images, labels, models, credentials or participant links are distributed in this repository. The data file is not supplied by this code license.

These exploratory sensitivity analyses were requested after examination of the main study. They are not new independent confirmation.
