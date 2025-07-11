import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) when mapping specific
    German "mouldy" concepts to English while maintaining discrimination.
    """
    # Step 1: Define the original feature sets for the German concepts.
    # We create plausible feature sets for each concept.
    german_mouldy_cellar_features = {'mould_present', 'dampness', 'stale_air', 'earthy_notes', 'cold_stone_sensation'}
    german_mouldy_fabric_features = {'mould_present', 'dampness', 'stale_air', 'decaying_textile_notes'}

    num_original_cellar_features = len(german_mouldy_cellar_features)
    num_original_fabric_features = len(german_mouldy_fabric_features)

    # Step 2: Determine the minimum number of features to retain.
    # To maintain the core concept of "mouldy", we must retain at least one common feature.
    # The common, or "core", features are:
    core_features = german_mouldy_cellar_features.intersection(german_mouldy_fabric_features)
    # Let's assume we must retain at least one of these, e.g., 'mould_present'.
    min_retained_core_features = 1

    # To distinguish between 'cellar' and 'fabric', we must also retain at least
    # one feature that is unique to each concept.
    min_retained_distinguishing_features = 1

    # Therefore, the minimum number of retained features for each concept is the sum.
    num_retained_features = min_retained_core_features + min_retained_distinguishing_features

    # Step 3: Calculate the FPR for each concept.
    # FPR = Retained Features / Original Features
    fpr_cellar = num_retained_features / num_original_cellar_features
    fpr_fabric = num_retained_features / num_original_fabric_features

    # Step 4: Find the minimum FPR achievable in this mapping scenario.
    min_fpr = min(fpr_cellar, fpr_fabric)

    # Explain the result step-by-step
    print("--- Conceptual Model ---")
    print(f"Original features for 'mouldy cellar': {num_original_cellar_features} {german_mouldy_cellar_features}")
    print(f"Original features for 'mouldy fabric': {num_original_fabric_features} {german_mouldy_fabric_features}")
    print("\n--- Minimum Retention for Discrimination ---")
    print(f"To preserve the core 'mouldy' concept, we must retain at least {min_retained_core_features} core feature (e.g., 'mould_present').")
    print(f"To distinguish between the two, we must retain at least {min_retained_distinguishing_features} unique feature for each.")
    print(f"Total minimum retained features per concept = {num_retained_features}")
    print("\n--- FPR Calculation ---")
    print(f"FPR for 'mouldy cellar' = Retained Features / Original Features")
    print(f"FPR for 'mouldy cellar' = {num_retained_features} / {num_original_cellar_features} = {fpr_cellar:.2f}")
    print(f"\nFPR for 'mouldy fabric' = Retained Features / Original Features")
    print(f"FPR for 'mouldy fabric' = {num_retained_features} / {num_original_fabric_features} = {fpr_fabric:.2f}")
    print("\n--- Conclusion ---")
    print(f"The minimum achievable FPR is the smaller of the two calculated rates.")
    print(f"Minimum FPR = min({fpr_cellar:.2f}, {fpr_fabric:.2f}) = {min_fpr:.2f}")

    # Final answer format
    print("\n<<<" + str(min_fpr) + ">>>")


calculate_min_fpr()