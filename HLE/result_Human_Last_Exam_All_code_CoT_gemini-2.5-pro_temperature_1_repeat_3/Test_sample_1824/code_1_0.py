def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping
    'mouldy' from German to English while preserving key distinctions.
    """
    # Step 1: Define the granular German feature sets.
    # We model the features for 'mouldy cellar' and 'mouldy fabric'.
    # These include common features ('damp', 'fungal', 'decay') and specific nuance features.
    features_mouldy_cellar = {'damp', 'fungal', 'decay', 'stale_air', 'mineral_earthy'}
    features_mouldy_fabric = {'damp', 'fungal', 'decay', 'fabric_substrate'}

    # The total pool of original features is the union of all features from the specific concepts.
    original_features = features_mouldy_cellar.union(features_mouldy_fabric)
    num_original_features = len(original_features)

    # Step 2: Determine the minimum features to retain for the English mapping.
    # To be considered 'mouldy' at all, some core features must be kept.
    retained_core = {'damp', 'fungal'}

    # To distinguish 'cellar' from 'fabric', we must keep their unique identifying features.
    # The feature unique to the cellar concept in our model is 'mineral_earthy'.
    # The feature unique to the fabric concept is 'fabric_substrate'.
    retained_discriminating = {'mineral_earthy', 'fabric_substrate'}

    # The minimum set of retained features includes both the core and discriminating features.
    retained_features = retained_core.union(retained_discriminating)
    num_retained_features = len(retained_features)

    # Step 3: Calculate the Feature Retention Rate (FPR).
    fpr = num_retained_features / num_original_features

    # Step 4: Print the results and the final equation.
    print(f"Total original features in German concept space: {sorted(list(original_features))}")
    print(f"Minimum retained features in English to preserve meaning and discrimination: {sorted(list(retained_features))}\n")

    print("Calculation of Minimum Feature Retention Rate (FPR):")
    print(f"FPR = (Number of Retained Features) / (Number of Original Features)")
    print(f"The final equation is: {num_retained_features} / {num_original_features} = {fpr}")

calculate_minimum_fpr()
