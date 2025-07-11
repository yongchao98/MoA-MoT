import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while maintaining discrimination between contexts.
    """
    # Step 1: Model the features for the nuanced German concept of 'mouldy'.
    # Let's assume the German concept has 5 distinct features that capture its nuances.
    # f1: core musty scent
    # f2: damp earth/stone notes (common in cellars)
    # f3: general organic decay notes
    # f4: heavy, stale air (common in cellars)
    # f5: sour/sharp notes from fabric decay
    german_original_features = {'f1_musty', 'f2_damp_earth', 'f3_decay', 'f4_heavy_stale_air', 'f5_fabric_sour'}
    num_original_features = len(german_original_features)

    # Step 2: Define which features apply to each specific context.
    mouldy_cellar_features = {'f1_musty', 'f2_damp_earth', 'f3_decay', 'f4_heavy_stale_air'}
    mouldy_fabric_features = {'f1_musty', 'f3_decay', 'f5_fabric_sour'}

    # Step 3: Identify core features (common to both contexts).
    # The mapping must preserve these to be a meaningful translation of "mouldy".
    core_features = mouldy_cellar_features.intersection(mouldy_fabric_features)
    num_core_features = len(core_features)

    # Step 4: Apply the discrimination constraint.
    # To discriminate between 'cellar' and 'fabric', the English concept must retain
    # at least one feature that is NOT common to both.
    # If only core features were retained, both 'mouldy cellar' and 'mouldy fabric'
    # would map to the exact same set of features, making them indistinguishable.
    # Therefore, we must retain all core features plus a minimum of one discriminating feature.
    min_retained_features = num_core_features + 1

    # Step 5: Calculate the minimum FPR.
    min_fpr = min_retained_features / num_original_features

    print("This script calculates the minimum Feature Retention Rate (FPR).")
    print("-" * 60)
    print(f"1. Total features in the original German concept: {num_original_features}")
    print(f"2. Core features common to both 'cellar' and 'fabric': {num_core_features}")
    print("3. To maintain discrimination, we must retain all core features plus at least one non-core feature.")
    print(f"   Minimum features to retain = {num_core_features} (core) + 1 (discriminating) = {min_retained_features}")
    print("-" * 60)
    print("Final Calculation:")
    print(f"Minimum FPR = (Minimum Retained Features) / (Original Features)")
    print(f"Minimum FPR = {min_retained_features} / {num_original_features} = {min_fpr}")

calculate_min_fpr()
<<<0.6>>>