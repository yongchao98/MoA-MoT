# 1. Model the original German semantic features.
# Let's assume the granular German concept 'moderig' (mouldy) is composed of:
# - A number of core features for the general "mouldy" concept.
# - A number of specific features that describe a "mouldy cellar".
# - A number of specific features that describe "mouldy fabric".

core_features = 5
specific_features_cellar = 2
specific_features_fabric = 3

# The total number of original features is the sum of all distinct features.
original_features_total = core_features + specific_features_cellar + specific_features_fabric

# 2. Determine the minimum features to retain for the English mapping.
# To preserve the core meaning of "mouldy", we must retain at least one core feature.
retained_core_min = 1

# To discriminate between "mouldy cellar" and "mouldy fabric", we must
# retain at least one specific feature from each category.
retained_cellar_min = 1
retained_fabric_min = 1

# The minimum total number of retained features is the sum of these minimums.
retained_features_total = retained_core_min + retained_cellar_min + retained_fabric_min

# 3. Calculate the Feature Retention Rate (FPR).
fpr = retained_features_total / original_features_total

# Output the final equation and the result.
print("To solve this, we define the feature counts:")
print(f"Original Core Features: {core_features}")
print(f"Original Specific 'Cellar' Features: {specific_features_cellar}")
print(f"Original Specific 'Fabric' Features: {specific_features_fabric}")
print("-" * 20)
print(f"Minimum Retained Core Features: {retained_core_min}")
print(f"Minimum Retained 'Cellar' Features: {retained_cellar_min}")
print(f"Minimum Retained 'Fabric' Features: {retained_fabric_min}")
print("-" * 20)
print("Final Calculation:")
print(f"Total Retained Features = {retained_core_min} + {retained_cellar_min} + {retained_fabric_min} = {retained_features_total}")
print(f"Total Original Features = {core_features} + {specific_features_cellar} + {specific_features_fabric} = {original_features_total}")
print(f"Minimum FPR = {retained_features_total} / {original_features_total} = {fpr}")
