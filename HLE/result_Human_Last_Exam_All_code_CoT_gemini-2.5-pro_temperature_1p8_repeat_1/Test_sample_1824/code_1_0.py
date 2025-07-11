import math

# Step 1: Model the concepts with hypothetical semantic features.
# We create plausible feature sets for the two specific German olfactory concepts.
german_mouldy_cellar_features = {'damp', 'earthy', 'stale_air', 'decaying_wood', 'confined_space', 'high_humidity'}
german_mouldy_fabric_features = {'damp', 'stale_air', 'decaying_fabric', 'chemical_undertone', 'high_humidity'}

print("Step 1: Define hypothetical feature sets for German concepts.")
print(f"Features for 'mouldy cellar': {german_mouldy_cellar_features}")
print(f"Features for 'mouldy fabric': {german_mouldy_fabric_features}")
print("-" * 20)

# Step 2: Define the total set of "Original Features".
# This is the union of all features from both specific concepts, representing the
# richness of the general German concept 'moderig' (mouldy).
original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
num_original_features = len(original_features)

print("Step 2: Calculate the total number of original features.")
print(f"The complete set of original features is the union of both sets: {original_features}")
print(f"Total number of Original Features = {num_original_features}")
print("-" * 20)

# Step 3: Determine the minimum features needed for discrimination.
# To discriminate, the mapped feature sets must not be identical.
# First, find the core features shared by both concepts.
core_features = german_mouldy_cellar_features.intersection(german_mouldy_fabric_features)
num_core_features = len(core_features)

print("Step 3: Determine the minimum features to retain for discrimination.")
print(f"Core features shared by both concepts: {core_features}")
print(f"Number of core features = {num_core_features}")
print("\nIf only these core features are retained, we cannot distinguish between the two smells.")
print("To maintain discrimination, we must retain all core features plus at least one distinguishing feature.")

# The minimum number of features to retain is the number of core features plus one.
min_retained_features = num_core_features + 1
print(f"Minimum Retained Features = (Number of core features) + 1 = {num_core_features} + 1 = {min_retained_features}")
print("-" * 20)

# Step 4: Calculate the minimum Feature Retention Rate (FPR).
fpr = min_retained_features / num_original_features

print("Step 4: Calculate the minimum Feature Retention Rate (FPR).")
print("FPR = Retained Features / Original Features")
print(f"FPR = {min_retained_features} / {num_original_features}")
print(f"The minimum FPR is: {fpr:.4f}")

# Final Answer in the required format
final_answer = min_retained_features / num_original_features
print(f"<<<{final_answer}>>>")