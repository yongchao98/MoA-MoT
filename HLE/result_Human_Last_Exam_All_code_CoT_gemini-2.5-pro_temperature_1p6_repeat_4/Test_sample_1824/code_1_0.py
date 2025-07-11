import math

# Step 1: Define the granular German feature sets.
# Let's assume 'mouldy cellar' has features for mould, dampness, earth, stale air, and stone.
# Let's assume 'mouldy fabric' has features for mould, textile, stale air, and storage chemicals.
german_mouldy_cellar_features = {'f_mould', 'f_damp', 'f_earthy', 'f_stale_air', 'f_stone'}
german_mouldy_fabric_features = {'f_mould', 'f_textile', 'f_stale_air', 'f_chemical'}

num_original_cellar = len(german_mouldy_cellar_features)
num_original_fabric = len(german_mouldy_fabric_features)

# Step 2 & 3: Model the minimal feature retention in English to maintain discrimination.
# To maintain discrimination with the fewest possible features, the two English concept sets
# must not be identical. The most efficient way to achieve this is to retain the core
# 'mouldy' feature in both, and add just one unique feature to one concept. This results
# in one concept being a subset of the other, which is the minimal possible difference.

# For example, let's keep {'f_mould'} for 'mouldy fabric'.
# To distinguish it, 'mouldy cellar' needs at least one more feature, e.g. {'f_mould', 'f_earthy'}.
# The sets {'f_mould'} and {'f_mould', 'f_earthy'} are distinct.

min_retained_fabric = 1 # Representing {'f_mould'}
min_retained_cellar = 2 # Representing {'f_mould', 'f_earthy'}

# Step 4: Calculate the Feature Retention Rate (FPR) for both concepts.
fpr_cellar = min_retained_cellar / num_original_cellar
fpr_fabric = min_retained_fabric / num_original_fabric

# Step 5: Determine the minimum achievable FPR.
# The minimum FPR achievable is the lesser of the two rates, as it's the lowest
# rate within the system that successfully maintains discrimination.
min_fpr = min(fpr_cellar, fpr_fabric)

# Identify which calculation corresponds to the minimum FPR.
if min_fpr == fpr_fabric:
    retained_features = min_retained_fabric
    original_features = num_original_fabric
    concept_name = "'mouldy fabric'"
else:
    retained_features = min_retained_cellar
    original_features = num_original_cellar
    concept_name = "'mouldy cellar'"

print(f"To distinguish 'mouldy cellar' from 'mouldy fabric' with minimum features:")
print(f" - Original features for 'mouldy cellar': {num_original_cellar}")
print(f" - Original features for 'mouldy fabric': {num_original_fabric}")
print(f" - Minimum retained features for 'mouldy fabric': {min_retained_fabric}")
print(f" - Minimum retained features for 'mouldy cellar': {min_retained_cellar}")
print("-" * 20)
print(f"FPR for 'mouldy fabric' = {min_retained_fabric} / {num_original_fabric} = {fpr_fabric}")
print(f"FPR for 'mouldy cellar' = {min_retained_cellar} / {num_original_cellar} = {fpr_cellar}")
print("-" * 20)
print("The minimum achievable FPR is the smaller of these two values.")
print("The final calculation for the minimum FPR is:")
print(f"{retained_features} / {original_features} = {min_fpr}")

# Final answer block as requested by the user format
# The value is min_fpr
# print(f'<<<{min_fpr}>>>')