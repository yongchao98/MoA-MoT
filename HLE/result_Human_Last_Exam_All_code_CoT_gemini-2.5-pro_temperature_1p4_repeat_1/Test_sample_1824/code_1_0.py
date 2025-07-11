def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) based on a model
    of semantic features for olfactory concepts.
    """

    # Step 1: Define the number of original features in the high-granularity language (German).
    # We model this based on the problem description:
    # - 1 core feature shared by all 'mouldy' concepts.
    # - 2 specific features for the 'mouldy cellar' smell (e.g., 'damp', 'earthy').
    # - 2 specific features for the 'mouldy fabric' smell (e.g., 'stale', 'dusty').
    # Total original features = 1 (core) + 2 (cellar) + 2 (fabric).
    original_features = 5

    # Step 2: Define the minimum number of features that must be retained in the
    # low-granularity language (English) to maintain discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric', we need:
    # - The 1 core 'mouldy' feature.
    # - At least 1 representative feature for the 'cellar' type.
    # - At least 1 representative feature for the 'fabric' type.
    # Total minimum retained features = 1 + 1 + 1.
    retained_features = 3

    # Step 3: Calculate the Feature Retention Rate (FPR).
    fpr = retained_features / original_features

    # Step 4: Output the explanation and the final equation with all numbers.
    print(f"Modeling the semantic features:")
    print(f"Total Original Features (German) = {original_features}")
    print(f"Minimum Retained Features (English) for discrimination = {retained_features}")
    print("\nThe minimum Feature Retention Rate (FPR) is Retained Features / Original Features.")
    print(f"Final Equation: FPR = {retained_features} / {original_features}")
    print(f"Result: {fpr}")

calculate_minimum_fpr()
<<<0.6>>>