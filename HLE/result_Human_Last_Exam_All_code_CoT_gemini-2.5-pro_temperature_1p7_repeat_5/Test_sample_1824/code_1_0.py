import sys

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """
    # Step 1: Define hypothetical feature sets for the German concepts.
    # The German word for 'mouldy' can be nuanced (e.g., 'moderig').
    # Let's model the features for specific contexts.
    german_mouldy_cellar = {'musty', 'damp', 'earthy', 'cold', 'stone-like'}
    german_mouldy_fabric = {'musty', 'damp', 'stale', 'dusty', 'organic-decay'}

    # Step 2: Calculate the total number of original features.
    # This is the union of all features used to describe the concepts in German.
    original_features_set = german_mouldy_cellar.union(german_mouldy_fabric)
    num_original_features = len(original_features_set)

    # Step 3: Determine the minimum features needed to retain discrimination in English.
    # In English, both are just 'mouldy'. They must share a core feature.
    # Let's assume the core 'mouldy' feature is 'musty'.
    
    # To maintain discrimination, their feature sets in English cannot be identical.
    # We want the smallest possible total number of retained features.
    # The most minimal way is to retain the core feature for both, and add one
    # distinguishing feature to just one of the concepts.

    # English concept for 'mouldy fabric' retains only the core feature.
    retained_fabric_features = {'musty'}
    
    # English concept for 'mouldy cellar' retains the core feature plus one
    # of its unique original features (e.g., 'earthy') to distinguish it.
    retained_cellar_features = {'musty', 'earthy'}

    # The total set of features retained in the English mapping is the union of these two minimal sets.
    total_retained_features_set = retained_fabric_features.union(retained_cellar_features)
    num_retained_features = len(total_retained_features_set)

    # Step 4: Calculate the minimum FPR.
    if num_original_features == 0:
        # Avoid division by zero, although it shouldn't happen with our sets.
        min_fpr = 0
    else:
        min_fpr = num_retained_features / num_original_features

    # Print the explanation and the final equation.
    print("This script models the semantic features of 'mouldy cellar' and 'mouldy fabric'.")
    print(f"Original features for 'mouldy cellar' in German: {german_mouldy_cellar}")
    print(f"Original features for 'mouldy fabric' in German: {german_mouldy_fabric}")
    print(f"Total unique original features: {original_features_set}")
    print(f"Number of Original Features = {num_original_features}\n")

    print("To maintain discrimination with minimal feature loss:")
    print(f" - Retained English features for 'mouldy fabric': {retained_fabric_features}")
    print(f" - Retained English features for 'mouldy cellar': {retained_cellar_features}")
    print(f"Total unique retained features for the mapping: {total_retained_features_set}")
    print(f"Minimum Number of Retained Features = {num_retained_features}\n")

    print("The minimum Feature Retention Rate (FPR) is calculated as:")
    print(f"FPR = Retained Features / Original Features")
    print(f"FPR = {num_retained_features} / {num_original_features} = {min_fpr:.2f}")
    
    # Final answer in the required format
    # Redirecting the final answer to a different stream to avoid it being mixed with print statements
    # The format requested is <<<answer>>>, so we will print that last.
    sys.stderr.write(f"<<<{min_fpr}>>>")

calculate_min_fpr()