import math

def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping a granular
    German olfactory concept to English while preserving specific distinctions.
    """

    # Step 1 & 2: Define feature sets for the olfactory concepts.
    # The German concept 'moderig' is broad and contains features for various contexts.
    german_moderig_features = {
        'f1_core_mould',      # The fundamental smell of mould
        'f2_dampness',        # The scent of moisture, common to all
        'f3_earthy',          # Soil, mineral-like scent (strong in cellars)
        'f4_stale_air',       # Lack of ventilation
        'f5_organic_decay',   # Rotting wood/leaves (associated with cellars)
        'f6_chemical_textile' # Odor associated with decaying fabric/dyes
    }

    # The English concept 'mouldy cellar' is a specific context.
    english_mouldy_cellar = {
        'f1_core_mould',
        'f2_dampness',
        'f3_earthy',
        'f4_stale_air',
        'f5_organic_decay'
    }

    # The English concept 'mouldy fabric' is another specific context.
    english_mouldy_fabric = {
        'f1_core_mould',
        'f2_dampness',
        'f4_stale_air',
        'f6_chemical_textile'
    }

    # Step 3: Identify core and distinguishing features.
    # We must retain the core concept of 'mould' to preserve meaning.
    core_features = {'f1_core_mould'}

    # To distinguish between 'cellar' and 'fabric', we need the features
    # that are not shared between them.
    distinguishing_features_cellar_only = english_mouldy_cellar.difference(english_mouldy_fabric)
    distinguishing_features_fabric_only = english_mouldy_fabric.difference(english_mouldy_cellar)
    
    # We need to pick the minimum number of features to preserve the distinction.
    # This means at least one from each "only" set.
    # Let's pick 'f3_earthy' for cellar and 'f6_chemical_textile' for fabric.
    min_distinguishing_set = {'f3_earthy', 'f6_chemical_textile'}


    # Step 4: Calculate the minimum set of features that must be retained.
    # This is the union of the core feature(s) and the minimal distinguishing features.
    min_retained_features = core_features.union(min_distinguishing_set)

    num_original_features = len(german_moderig_features)
    num_retained_features = len(min_retained_features)

    # Step 5: Calculate the minimum Feature Retention Rate (FPR).
    min_fpr = num_retained_features / num_original_features

    print("--- Analysis ---")
    print(f"Original German features ('moderig'): {german_moderig_features}")
    print(f"Total original features: {num_original_features}\n")
    
    print(f"Features for 'mouldy cellar': {english_mouldy_cellar}")
    print(f"Features for 'mouldy fabric': {english_mouldy_fabric}\n")

    print(f"Features unique to 'cellar': {distinguishing_features_cellar_only}")
    print(f"Features unique to 'fabric': {distinguishing_features_fabric_only}\n")

    print(f"To preserve the basic meaning, we must retain: {core_features}")
    print(f"To distinguish 'cellar' from 'fabric', we must additionally retain at least: {min_distinguishing_set}")
    print(f"Minimum set of features to retain: {min_retained_features}")
    print(f"Total minimum retained features: {num_retained_features}\n")

    print("--- Calculation ---")
    print("FPR = (Minimum Retained Features) / (Original Features)")
    # The final equation with the numbers plugged in
    print(f"FPR = {num_retained_features} / {num_original_features} = {min_fpr:.1f}")

    return min_fpr

# Execute the function and capture the final answer.
final_answer = calculate_minimum_fpr()
print(f"\n<<<The minimum Feature Retention Rate (FPR) is {final_answer:.1f}>>>")
print(f"<<<{final_answer:.1f}>>>")