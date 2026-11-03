# Reduced-data schema

Top level: groups, lambda_original, lambda_cell_reference, bootstrap (seed and draws).
Each group: group_key (opaque local key), arms, estimated_baseline_dice (WT,TC,ET,NCR/NETC), candidates.
Each arm A,B,C,point_B: counts (four [predicted volume, reference volume, intersection] integer triplets); transitions (changed, repaired, new_errors, wrong_to_wrong).
Each candidate: mass, source, target, alpha (four positive concentrations), gt_tp_increment (four signed integers), transitions (repaired,new_errors,wrong_to_wrong).
Computational labels: background0, NCR/NETC1,SNFH2,ET3; WT={1,2,3},TC={1,3},ET={3},NCR/NETC={1}.
All 256 candidates per group refer to disjoint fixed cells. No image coordinate, acquisition identifier or original file hash is needed for these numerical checks. Opaque keys are pseudonyms; their absence from a public package is deliberate.
