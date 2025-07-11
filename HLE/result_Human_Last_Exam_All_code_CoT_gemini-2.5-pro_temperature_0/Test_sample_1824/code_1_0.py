import math

# Step 1: Define the granular set of features for the original German concept.
# This represents the rich understanding of 'mouldy' in German.
original_german_features = {
    # Core features essential for any 'mouldy' concept
    'f1_has_fungus', 'f2_is_damp', 'f3_is_earthy',
    # Features specific to 'mouldy cellar'
    'f4_has_mineral_notes', 'f5_has_stagnant_air',
    # Features specific to 'mouldy fabric'
    'f6_has_decaying_textile', 'f7_is_in_confined_space'
}
num_original_features = len(original_german_features)

# Step 2: Identify the minimum set of features to retain for the English concept.
# To preserve the core meaning, we must keep the core features.
retained_core_features = {'f1_has_fungus', 'f2_is_damp', 'f3_is_earthy'}

# To maintain discrimination, we must keep at least one unique feature for each sub-type.
# We select the most representative feature for each.
retained_discriminating_features = {'f4_has_mineral_notes', 'f6_has_decaying_textile'}

# The minimum set of retained features is the union of these two sets.
minimum_retained_features = retained_core_features.union(retained_discriminating_features)
num_retained_features = len(minimum_retained_features)

# Step 3: Calculate the Feature Retention Rate (FPR).
fpr = num_retained_features / num_original_features

# Output the final equation and the result.
print("Problem Analysis:")
print(f"1. Original German concept has {num_original_features} features: {original_german_features}")
print(f"2. To preserve the core meaning and discrimination, a minimum of {num_retained_features} features must be retained: {minimum_retained_features}")
print("\nCalculation:")
print(f"Minimum Feature Retention Rate (FPR) = (Number of Retained Features) / (Number of Original Features)")
# The final print statement shows the numbers used in the equation as requested.
print(f"Minimum FPR = {num_retained_features} / {num_original_features} = {fpr:.4f}")
