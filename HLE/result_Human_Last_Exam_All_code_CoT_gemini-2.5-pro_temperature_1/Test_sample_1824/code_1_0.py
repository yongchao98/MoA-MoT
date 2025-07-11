import sys

def solve_fpr_problem():
    """
    This script calculates the minimum Feature Retention Rate (FPR) for mapping
    a German 'mouldy' concept to English while preserving discrimination between
    'mouldy cellar' and 'mouldy fabric'.
    """

    # Step 1: Define the original, high-granularity semantic features for the German concepts.
    # We hypothesize a plausible set of features based on the concepts.
    features_cellar_de = {'mould', 'damp', 'earthy', 'stale_air', 'cold'}
    features_fabric_de = {'mould', 'damp', 'stale_air', 'fabric_smell'}

    print("--- Step 1: Defining Original Semantic Features ---")
    print(f"Features for 'mouldy cellar': {features_cellar_de}")
    print(f"Features for 'mouldy fabric': {features_fabric_de}\n")

    # Step 2: Determine the total set and number of original features.
    # This is the union of all features used across the concepts.
    total_original_features = features_cellar_de.union(features_fabric_de)
    num_original_features = len(total_original_features)

    print("--- Step 2: Calculating Total Original Features ---")
    print(f"The complete set of original features is: {total_original_features}")
    print(f"Total number of Original Features = {num_original_features}\n")

    # Step 3: Determine the minimum number of features to retain.
    # Constraint 1: The core concept 'mouldy' must be preserved, so 'mould' is retained.
    # Constraint 2: 'mouldy cellar' and 'mouldy fabric' must be distinguishable.
    # To satisfy this with the minimum number of features, we must retain the core 'mould' feature
    # plus at least one feature that differentiates the two concepts (e.g., 'earthy', which is
    # only in the cellar concept).
    # Minimal retained features for cellar: {'mould', 'earthy'}
    # Minimal retained features for fabric: {'mould'}
    # The two sets are now distinct. The total pool of features we needed to retain is {'mould', 'earthy'}.
    min_retained_features_set = {'mould', 'earthy'}
    num_retained_features = len(min_retained_features_set)
    
    print("--- Step 3: Determining Minimum Retained Features for Discrimination ---")
    print("To distinguish the concepts, their feature sets in English must differ.")
    print("The most efficient way is to retain the core feature ('mould') and one unique feature (e.g., 'earthy').")
    print(f"Minimum set of features to retain in the system: {min_retained_features_set}")
    print(f"Minimum number of Retained Features = {num_retained_features}\n")

    # Step 4: Calculate the minimum Feature Retention Rate (FPR).
    if num_original_features == 0:
        min_fpr = 0.0
    else:
        min_fpr = num_retained_features / num_original_features

    print("--- Step 4: Calculating the Minimum FPR ---")
    print("FPR = Retained Features / Original Features")
    print(f"FPR = {num_retained_features} / {num_original_features}")
    print(f"Minimum FPR = {min_fpr:.4f}\n")
    
    # Final answer for automated extraction
    # The format is specified by the user prompt.
    sys.stdout.write(f"<<<{min_fpr:.4f}>>>")

solve_fpr_problem()