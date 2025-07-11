def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping a German
    olfactory concept to English while preserving key distinctions.
    """
    # Step 1: Define the hypothetical set of original features for the German
    # concept 'moderig' (mouldy).
    original_german_features = {'damp', 'earthy', 'organic_decay', 'spores'}
    num_original_features = len(original_german_features)

    # Step 2: Identify the minimum features needed to discriminate between
    # 'mouldy cellar' (relies on 'earthy') and 'mouldy fabric' (relies on 'organic_decay').
    # To tell them apart, both of these discriminating features must be retained.
    retained_features_for_discrimination = {'earthy', 'organic_decay'}
    num_retained_features = len(retained_features_for_discrimination)

    # Step 3: Calculate the Feature Retention Rate (FPR).
    fpr = num_retained_features / num_original_features

    # Step 4: Print the results, showing the equation.
    print("Problem Analysis:")
    print(f"1. Original German features for 'mouldy': {original_german_features}")
    print(f"2. Features needed to distinguish 'cellar' vs. 'fabric': {retained_features_for_discrimination}")
    print("-" * 20)
    print("Calculation:")
    print(f"Total number of original features = {num_original_features}")
    print(f"Minimum number of retained features = {num_retained_features}")
    print("\nFinal Equation:")
    print(f"Minimum FPR = (Minimum Retained Features) / (Original Features)")
    print(f"Minimum FPR = {num_retained_features} / {num_original_features} = {fpr}")

calculate_minimum_fpr()
<<<0.5>>>