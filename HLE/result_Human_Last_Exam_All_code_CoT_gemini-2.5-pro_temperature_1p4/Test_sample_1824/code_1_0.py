import math

# Step 1: Define the original features based on German's high granularity.
# The total German conceptual space for mustiness is modeled as the union of features
# from its specific terms: 'moderig' (mouldy), 'muffig' (fusty), and 'stickig' (stuffy).
moderig_features = {'is_decaying', 'is_damp'}
muffig_features = {'is_stale', 'is_on_soft_material'}
stickig_features = {'is_in_enclosed_space', 'has_mineral_note'}

# The total set of original features is the union of all these nuanced features.
original_features = moderig_features.union(muffig_features).union(stickig_features)
num_original_features = len(original_features)

# Step 2: Model the specific German concepts by selecting from the original pool.
# 'mouldy cellar' combines mould ('moderig') and stuffiness ('stickig').
german_mouldy_cellar = {'is_decaying', 'is_damp', 'is_in_enclosed_space', 'has_mineral_note'}

# 'mouldy fabric' combines mould ('moderig') and fustiness ('muffig').
german_mouldy_fabric = {'is_decaying', 'is_damp', 'is_stale', 'is_on_soft_material'}

# Step 3: Determine the minimum features to retain for English 'mouldy'.
# The English concept must retain enough features to preserve the core meaning of "mouldy"
# AND the ability to distinguish the two contexts.
# The common "mouldy" features are {'is_decaying', 'is_damp'}.
# The features that discriminate between cellar and fabric are {'is_in_enclosed_space', 'has_mineral_note'} vs {'is_stale', 'is_on_soft_material'}.
# To maintain discrimination, the retained set must make the two final concepts different and non-empty.
# A) Retaining just 1 feature is not enough. For example, retaining only a common feature like 'is_decaying' would make both contexts map to {'is_decaying'}, losing discrimination. Retaining only a unique feature like 'has_mineral_note' would make the fabric context map to an empty set, losing its meaning.
# B) Therefore, we must retain at least 2 features: one common feature to preserve the core "mouldy" sense in both contexts, and one discriminating feature to tell them apart.
# For example, let's retain {'is_decaying', 'has_mineral_note'}.
# - Mapped cellar: {'is_decaying', 'has_mineral_note'}.intersection(german_mouldy_cellar) -> {'is_decaying', 'has_mineral_note'}
# - Mapped fabric: {'is_decaying', 'has_mineral_note'}.intersection(german_mouldy_fabric) -> {'is_decaying'}
# The resulting concepts are {'is_decaying', 'has_mineral_note'} and {'is_decaying'}. They are distinct and non-empty.
# Thus, the minimum number of retained features is 2.

num_retained_features = 2

# Step 4: Calculate and print the Feature Retention Rate (FPR).
feature_retention_rate = num_retained_features / num_original_features

print("--- Calculation Steps ---")
print(f"1. Total Original Features from German granular concepts ('moderig', 'muffig', 'stickig'): {num_original_features}")
print(f"2. Minimum Retained Features in English to maintain discrimination: {num_retained_features}")
print("\n--- Final Equation ---")
print(f"Minimum FPR = (Minimum Retained Features) / (Total Original Features)")
print(f"Minimum FPR = {num_retained_features} / {num_original_features}")
print(f"Result = {feature_retention_rate}")
print(f"<<<{feature_retention_rate}>>>")