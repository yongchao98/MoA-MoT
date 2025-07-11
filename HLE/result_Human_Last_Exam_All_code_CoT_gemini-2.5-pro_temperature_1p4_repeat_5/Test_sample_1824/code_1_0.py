def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """

    # Step 1: Model the rich German concept 'moderig' (mouldy) with a set of features.
    # German has high granularity, so we define a comprehensive set of 5 features.
    original_german_features = {
        'f1_mould_presence',  # The core concept: something is mouldy.
        'f2_dampness',        # A feature of high humidity.
        'f3_stale_air',       # A feature for lack of ventilation.
        'f4_earthy_mineral',  # A feature specific to cellars, soil, stone.
        'f5_decaying_textile' # A feature specific to mould on fabric.
    }
    num_original_features = len(original_german_features)

    # Step 2: Represent the specific contexts as subsets of the original features.
    # 'mouldy cellar' includes the earthy/mineral note but not the textile note.
    mouldy_cellar_features = {'f1_mould_presence', 'f2_dampness', 'f3_stale_air', 'f4_earthy_mineral'}
    # 'mouldy fabric' includes the textile note but not the earthy/mineral note.
    mouldy_fabric_features = {'f1_mould_presence', 'f2_dampness', 'f3_stale_air', 'f5_decaying_textile'}

    # Step 3: Determine the minimum set of features to retain for discrimination.
    # The English mapping to 'mouldy' must retain the core feature 'f1_mould_presence'.
    # If only f1 is retained, both 'cellar' and 'fabric' would map to the same single feature,
    # and discrimination would be lost.
    # To distinguish them, we must also retain at least one of the features that make them different.
    # The differentiating features are 'f4_earthy_mineral' and 'f5_decaying_textile'.
    # Therefore, the minimum set of retained features must include the core feature plus one
    # of the differentiating features.
    # e.g., retained_set = {'f1_mould_presence', 'f4_earthy_mineral'}
    # With this set, 'cellar' maps to {f1, f4} and 'fabric' maps to {f1}. They are now distinct.
    min_retained_features = 2

    # Step 4: Calculate the Feature Retention Rate (FPR).
    fpr = min_retained_features / num_original_features

    # Step 5: Output the logic and the final calculation.
    print("--- Problem Analysis ---")
    print(f"1. Total number of original features in the German concept: {num_original_features}")
    print(f"2. Minimum features to retain for discrimination (core concept + 1 distinguishing feature): {min_retained_features}")
    print("\n--- Calculation ---")
    print("Feature Retention Rate (FPR) = Minimum Retained Features / Original Features")
    print(f"FPR = {min_retained_features} / {num_original_features}")
    print(f"Final FPR = {fpr}")

calculate_minimum_fpr()
<<<0.4>>>