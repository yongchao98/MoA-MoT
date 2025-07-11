def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping a granular
    'mouldy' concept while preserving discrimination between contexts.
    """

    # Step 1: Model the original, granular feature set for the German concept.
    # Let's assume the German concept of 'mouldy' ('moderig') has features that
    # cover the core smell as well as different contextual notes.
    # Features = {core_mould, damp, decay, earthy, mineral, textile}
    # A 'mouldy cellar' would activate {core_mould, damp, decay, earthy, mineral}.
    # A 'mouldy fabric' would activate {core_mould, damp, decay, textile}.
    original_features_count = 6

    # Step 2 & 3: Determine the minimum set of features to retain for discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric', the English mapping must
    # preserve the features that make them different, plus a core feature to
    # identify the concept as "mouldy".

    # - Retain one core feature to identify the concept: e.g., 'core_mould'. (1 feature)
    # - Retain the feature that distinguishes a cellar's smell: e.g., 'earthy'. (1 feature)
    # - Retain the feature that distinguishes a fabric's smell: e.g., 'textile'. (1 feature)
    # The minimum set of retained features would be {core_mould, earthy, textile}.
    retained_features_count = 3

    # Step 4: Calculate the Feature Retention Rate (FPR).
    fpr = retained_features_count / original_features_count

    print("Calculating the Minimum Feature Retention Rate (FPR)...")
    print("-" * 50)
    print(f"1. Assumed number of features in the original granular German concept: {original_features_count}")
    print(f"2. Minimum features needed to retain for discriminating 'mouldy cellar' vs. 'mouldy fabric': {retained_features_count}")
    print("   (These are: one core 'mouldy' feature, one 'cellar' specific feature, and one 'fabric' specific feature)")
    print("-" * 50)
    print("The final calculation is:")
    print(f"FPR = Retained Features / Original Features")
    print(f"FPR = {retained_features_count} / {original_features_count} = {fpr}")


if __name__ == "__main__":
    calculate_minimum_fpr()
    print("\n<<<0.5>>>")
