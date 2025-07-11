import math

# Step 1: Model the concepts with semantic features.
# Let's assume the core concept 'mouldy' has 4 features.
core_mouldy_features = {'f_damp', 'f_decay', 'f_earthy', 'f_fungal'}

# Let's assume the context 'cellar' adds 2 specific nuance features.
cellar_nuance_features = {'f_mineral', 'f_stagnant_air'}

# And the context 'fabric' adds 2 different nuance features.
fabric_nuance_features = {'f_textile', 'f_close_proximity'}

# The full German concept for 'mouldy cellar' combines core and nuance features.
german_mouldy_cellar = core_mouldy_features.union(cellar_nuance_features)

# The full German concept for 'mouldy fabric' combines core and nuance features.
german_mouldy_fabric = core_mouldy_features.union(fabric_nuance_features)

# Step 2: Determine the total number of original features.
# This is the union of all features required to describe both distinct concepts.
total_original_features_set = german_mouldy_cellar.union(german_mouldy_fabric)
num_original_features = len(total_original_features_set)

# Step 3: Determine the minimum number of features that must be retained.
# The problem states core features must be preserved.
num_retained_core = len(core_mouldy_features)

# To maintain discrimination, the mapped feature sets for cellar and fabric cannot be identical.
# If we only keep the core features, both concepts would map to the same set, losing discrimination.
# Therefore, we must retain at least one of the nuance features. The minimum is 1.
num_retained_nuance = 1

# The minimum number of retained features is the sum of core features and the minimum required nuance features.
min_num_retained_features = num_retained_core + num_retained_nuance

# Step 4: Calculate the minimum Feature Retention Rate (FPR).
# FPR = Retained Features / Original Features
min_fpr = min_num_retained_features / num_original_features

# Output the reasoning and the final calculation, showing each number.
print("To solve this, we define a model for the concepts:")
print(f"- Core 'mouldy' features: {num_retained_core}")
print(f"- 'Cellar' nuance features: {len(cellar_nuance_features)}")
print(f"- 'Fabric' nuance features: {len(fabric_nuance_features)}\n")

print(f"1. Calculate Total Original Features:")
print(f"   The total unique features for both concepts is {num_original_features}.\n")

print(f"2. Calculate Minimum Retained Features:")
print(f"   To preserve the core concept, we must retain all {num_retained_core} core features.")
print(f"   To distinguish 'cellar' from 'fabric', we must retain at least {num_retained_nuance} nuance feature.")
print(f"   Minimum features to retain = {num_retained_core} (core) + {num_retained_nuance} (nuance) = {min_num_retained_features}.\n")

print("3. Final Calculation of Minimum FPR:")
print(f"   Minimum FPR = Minimum Retained Features / Total Original Features")
print(f"   Minimum FPR = {min_num_retained_features} / {num_original_features} = {min_fpr}")
<<<0.625>>>