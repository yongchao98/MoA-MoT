def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between contexts.
    """
    # Step 1: Model the granular German concepts with semantic feature sets.
    # German has specific words for musty odors ('muffig', 'stickig', 'moderig').
    # We'll assign features based on these nuances.
    german_mouldy_cellar_features = {
        'has_mould',              # Core 'mouldy' concept
        'is_damp',                # Core 'mouldy' concept (muffig)
        'is_decomposing',         # Core 'mouldy' concept (moderig)
        'is_enclosed_space',      # Specific to cellar (stickig/stuffy)
        'has_earthy_smell'        # Specific to cellar
    }

    german_mouldy_fabric_features = {
        'has_mould',              # Core 'mouldy' concept
        'is_damp',                # Core 'mouldy' concept (muffig)
        'is_decomposing',         # Core 'mouldy' concept (moderig)
        'has_textile_smell'       # Specific to fabric
    }

    # Step 2: Identify the total number of original features (the union of all features).
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    # Step 3: Identify the core features that must be retained (the intersection).
    core_features = german_mouldy_cellar_features.intersection(german_mouldy_fabric_features)
    num_core_features = len(core_features)

    # Step 4 & 5: To discriminate, we must retain the core features plus at least one
    # distinguishing feature. The minimum addition is 1.
    # If we only retained core features, the English concepts would be identical.
    min_additional_features_for_discrimination = 1
    num_retained_features = num_core_features + min_additional_features_for_discrimination

    # Step 6: Calculate the Feature Retention Rate (FPR).
    fpr = num_retained_features / num_original_features

    print("--- Analysis of Feature Retention ---")
    print(f"1. Total unique features in German concepts (Original Features): {num_original_features}")
    print(f"2. Core features common to both concepts: {num_core_features}")
    print("3. To distinguish the two concepts, we must retain the core features plus at least one unique feature.")
    print(f"4. Minimum features to retain in English (Retained Features): {num_core_features} (core) + {min_additional_features_for_discrimination} (distinguishing) = {num_retained_features}")
    print("\n--- Final Calculation ---")
    print(f"Minimum FPR = Retained Features / Original Features")
    print(f"Minimum FPR = {num_retained_features} / {num_original_features} = {fpr:.3f}")

calculate_minimum_fpr()
<<<0.667>>>