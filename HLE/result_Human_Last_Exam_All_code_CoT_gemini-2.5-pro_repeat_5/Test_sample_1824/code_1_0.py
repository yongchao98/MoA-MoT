import math

# Step 1: Model the semantic features for the German concepts.
# These are hypothetical features for demonstration based on the problem description.
features_mouldy_cellar = {'musty', 'damp', 'decay', 'earthy', 'stagnant_air'}
features_mouldy_fabric = {'musty', 'damp', 'decay', 'fuzzy_texture'}

# Step 2: Identify the total pool of original features and the core features.
# The total original features are all unique features for both concepts.
original_features_total = features_mouldy_cellar.union(features_mouldy_fabric)
num_original_features = len(original_features_total)

# The core features are those shared by both concepts.
# The problem states core features are preserved in the mapping.
core_features = features_mouldy_cellar.intersection(features_mouldy_fabric)

# Step 3: Determine the minimum set of features that must be retained.
# To maintain discrimination, the feature sets for the two concepts after mapping
# must not be identical.
# If we only retain the core features, both concepts would map to the same set,
# and discrimination would be lost.
# Therefore, we must retain the core features PLUS at least one feature that
# is unique to one of the concepts. The minimum is to add just ONE such feature.

# Let's find the features that are not core (i.e., discriminating features).
discriminating_features = original_features_total.difference(core_features)

# To create the minimum retained set, we take the core features and add one
# discriminating feature.
min_retained_features_set = core_features.union({list(discriminating_features)[0]})
num_retained_features = len(min_retained_features_set)

# Step 4: Calculate the Feature Retention Rate (FPR).
fpr = num_retained_features / num_original_features

# --- Output the results ---
print(f"Assumed Features for 'Mouldy Cellar': {features_mouldy_cellar}")
print(f"Assumed Features for 'Mouldy Fabric': {features_mouldy_fabric}")
print("-" * 20)
print(f"Total Original Features in German ({num_original_features}): {original_features_total}")
print(f"Core Features to be Preserved: {core_features}")
print(f"Minimum Retained Features for Discrimination ({num_retained_features}): {min_retained_features_set}")
print("-" * 20)
print("To calculate the minimum Feature Retention Rate (FPR), we use the formula:")
print("FPR = (Minimum Retained Features) / (Total Original Features)")
print("\nThe final equation is:")
print(f"{num_retained_features} / {num_original_features}")
print(f"\nResult: {fpr}")

# The problem asks for the final answer in a specific format.
# The value is 4/6 = 0.666...
# We'll use this value for the final answer block.
final_answer = 2/3