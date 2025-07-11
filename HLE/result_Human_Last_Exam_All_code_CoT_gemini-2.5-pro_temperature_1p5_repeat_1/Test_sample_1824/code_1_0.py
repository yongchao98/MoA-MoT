import math

def solve_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """
    # Step 1: Define hypothetical semantic feature sets.
    # Based on the German examples, we can infer distinct features for different contexts.
    # 'muffig' (musty/stuffy) for a cellar implies dampness, poor ventilation, and an earthy smell.
    # 'moderig' (mouldy/decaying) for fabric implies dampness, organic decay, and fungal growth on a specific material.
    features_mouldy_cellar = {'dampness', 'stale_air', 'earthy_smell', 'fungal_growth', 'stone_substrate'}
    features_mouldy_fabric = {'dampness', 'organic_decay', 'sour_note', 'fungal_growth', 'textile_substrate'}

    # Step 2: Determine the total number of original features (the union of both sets).
    original_features = features_mouldy_cellar.union(features_mouldy_fabric)
    num_original_features = len(original_features)

    # Step 3: Identify the core features that must be preserved (the intersection).
    # "Cross-language mapping preserves core semantic features."
    core_features = features_mouldy_cellar.intersection(features_mouldy_fabric)
    num_core_features = len(core_features)

    # Step 4: Determine the minimum number of retained features to allow discrimination.
    # The retained set must include all core features.
    # If we only retain core features, both 'mouldy cellar' and 'mouldy fabric' would
    # map to the same feature set {'dampness', 'fungal_growth'}, losing discrimination.
    # Therefore, we must retain at least one additional feature that is unique to one
    # of the concepts to tell them apart.
    # Minimum Retained Features = Number of Core Features + 1 discriminating feature.
    num_retained_features = num_core_features + 1

    # Step 5: Calculate the minimum Feature Retention Rate (FPR).
    min_fpr = num_retained_features / num_original_features

    # --- Output the process and result ---
    print("A plausible breakdown of the problem is as follows:")
    print(f"\n1. Features for 'mouldy cellar' could be: {features_mouldy_cellar}")
    print(f"2. Features for 'mouldy fabric' could be: {features_mouldy_fabric}")

    print(f"\nFrom this, we determine the Original, Core, and Retained features:")
    print(f"- Total Original Features (Union of sets): {num_original_features}")
    print(f"- Core Features to be Preserved (Intersection of sets): {num_core_features}")
    print(f"- Minimum Retained Features (Core + 1 discriminating feature): {num_retained_features}")

    print("\nTo calculate the minimum Feature Retention Rate (FPR), we use the formula:")
    print("FPR = (Minimum Retained Features) / (Total Original Features)")
    print(f"\nFinal Equation: {num_retained_features} / {num_original_features} = {min_fpr}")

solve_fpr()
<<<0.375>>>