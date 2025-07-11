import sys

def solve_fpr_problem():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping a German
    olfactory concept to English while maintaining discrimination between two contexts.
    """
    # Step 1: Define feature sets for the two contexts.
    # These represent the concepts with full detail.
    features_mouldy_cellar = {'f_mould', 'f_damp', 'f_earth', 'f_stone'}
    features_mouldy_fabric = {'f_mould', 'f_damp', 'f_textile', 'f_decay'}

    print("Step 1: Defining the semantic features for each context.")
    print(f"Features for 'mouldy cellar': {features_mouldy_cellar}")
    print(f"Features for 'mouldy fabric': {features_mouldy_fabric}\n")

    # Step 2: Model the high-granularity German concept ('moderig') as the union of features.
    # This represents the total feature pool before any are lost in translation.
    original_german_features = features_mouldy_cellar.union(features_mouldy_fabric)
    num_original_features = len(original_german_features)

    print("Step 2: Modeling the original German concept ('moderig').")
    print("This is the union of all features from both contexts due to high lexical granularity.")
    print(f"Original German Features: {original_german_features}")
    print(f"Total number of Original Features: {num_original_features}\n")

    # Step 3: Identify core features that must be preserved.
    # These are the features common to both contexts.
    core_features = features_mouldy_cellar.intersection(features_mouldy_fabric)
    num_core_features = len(core_features)

    print("Step 3: Identifying the core semantic features.")
    print("Core features are common to both contexts and must be retained in the mapping.")
    print(f"Core Features: {core_features}")
    print(f"Number of Core Features: {num_core_features}\n")

    # Step 4: Determine the minimum number of features to retain for discrimination.
    # If we only retain core features, both contexts map to the same set, losing discrimination.
    # We must retain the core features plus at least one additional feature that is unique
    # to one of the contexts.
    min_retained_features = num_core_features + 1

    print("Step 4: Determining the minimum features needed to maintain discrimination.")
    print("If only the core features were retained, 'mouldy cellar' and 'mouldy fabric' would become identical, and discrimination would be lost.")
    print("Therefore, we must retain all core features plus at least one additional discriminating feature.")
    print(f"Minimum Retained Features = (Number of Core Features) + 1")
    print(f"Minimum Retained Features = {num_core_features} + 1 = {min_retained_features}\n")

    # Step 5: Calculate the minimum Feature Retention Rate (FPR).
    min_fpr = min_retained_features / num_original_features

    print("Step 5: Calculating the Minimum Feature Retention Rate (FPR).")
    print("FPR = Retained Features / Original Features")
    print(f"FPR = {min_retained_features} / {num_original_features}")
    print(f"Minimum FPR = {min_fpr}")

    # Final Answer Formatting
    # Redirect print output to stderr to hide it from the final answer block
    original_stdout = sys.stdout
    sys.stdout = sys.stderr
    print("\n--- End of reasoning ---")
    sys.stdout = original_stdout
    
    print(f"<<<{min_fpr}>>>")


solve_fpr_problem()