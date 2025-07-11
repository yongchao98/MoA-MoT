import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'moderig'
    to 'mouldy' while preserving discrimination between two contexts.
    """
    # Step 1: Model the semantic features of the granular German concept 'moderig'.
    # We hypothesize that 'moderig' contains a core feature for "mould" and
    # at least two nuance features for different contexts.
    features_moderig_german = {
        "f_core_mould",       # The essential feature for any "mouldy" concept.
        "f_nuance_cellar",    # A feature for the damp, earthy aspect often found in cellars.
        "f_nuance_fabric"     # A feature for the decaying textile aspect.
    }
    original_features_count = len(features_moderig_german)

    print("Step 1: Define the Original Features for the German concept 'moderig'.")
    print(f"The granular concept 'moderig' is modeled with {original_features_count} features: {features_moderig_german}")
    print(f"Therefore, Original Features = {original_features_count}\n")

    # Step 2: Determine the minimum retained features for discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric', their resulting
    # feature sets in English must be different.

    # If we only retain the core feature {'f_core_mould'}, both 'mouldy cellar'
    # and 'mouldy fabric' would map to the same set. Discrimination would fail.
    # To succeed, we must retain the core feature PLUS at least one nuance feature.
    min_retained_features_count = 2 # e.g., {"f_core_mould", "f_nuance_cellar"}

    print("Step 2: Determine the Minimum Retained Features for discrimination.")
    print("To distinguish between 'mouldy cellar' and 'mouldy fabric', we must retain enough features")
    print("so their resulting representations are not identical.")
    print("Retaining only the core 'mould' feature is not enough.")
    print("We must retain the core feature and at least one of the nuance features.")
    print(f"Therefore, Minimum Retained Features = {min_retained_features_count}\n")

    # Step 3: Calculate the minimum Feature Retention Rate (FPR).
    min_fpr = min_retained_features_count / original_features_count

    print("Step 3: Calculate the minimum Feature Retention Rate (FPR).")
    print("FPR = Retained Features / Original Features")
    print(f"Minimum FPR = {min_retained_features_count} / {original_features_count}")
    print(f"Result = {min_fpr:.3f}")
    
    # Final answer in the specified format
    # The value is 2/3.
    final_answer = 2/3
    print(f"\n<<<{'%.3f' % final_answer}>>>")

calculate_min_fpr()