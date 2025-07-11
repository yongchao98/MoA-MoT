import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) to discriminate
    between 'mouldy cellar' and 'mouldy fabric'.
    """

    # Step 1: Model the high-granularity concepts with semantic features.
    # Based on the prompt, German has high granularity. We define plausible feature sets.
    # 'moderig' (decay), 'stickig' (stagnant air), etc.
    german_mouldy_cellar_features = {'is_mouldy', 'is_damp', 'has_earthy_notes', 'has_stone_notes', 'has_stagnant_air'}
    german_mouldy_fabric_features = {'is_mouldy', 'is_damp', 'has_fabric_notes', 'has_organic_decay_notes'}

    # Step 2: Determine the total number of original features for the calculation.
    # This is the union of all features for the concepts we need to discriminate between.
    total_original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(total_original_features)

    # Step 3: Find the minimum number of features to retain for discrimination.
    # To map the concept 'mouldy' at all, the core feature 'is_mouldy' must be retained.
    # This is our baseline.
    num_retained_features = 1 # for 'is_mouldy'

    # Now, check if this is enough to discriminate.
    # If we only keep {'is_mouldy'}, both concepts would map to {'is_mouldy'}, and we cannot tell them apart.
    
    # To discriminate, we must retain at least one more feature that is unique to one of the concepts.
    # By adding one unique feature (e.g., 'has_stone_notes' or 'has_fabric_notes') to our retained set,
    # the mapped feature sets will become different.
    # e.g., mapping rule retains {'is_mouldy', 'has_stone_notes'}
    # Mapped cellar -> {'is_mouldy', 'has_stone_notes'}
    # Mapped fabric -> {'is_mouldy'}
    # These are now discriminable.
    
    # So, the minimum number of retained features is the core feature (1) plus one differentiating feature (1).
    num_retained_features += 1
    
    # Step 4: Calculate the Feature Retention Rate (FPR).
    fpr = num_retained_features / num_original_features

    print("--- Calculation of Minimum Feature Retention Rate (FPR) ---")
    print(f"Total set of original features: {total_original_features}")
    print(f"Number of Original Features = {num_original_features}")
    print(f"Minimum Retained Features to preserve meaning and allow discrimination = {num_retained_features}")
    print(f"The final equation is: FPR = Retained Features / Original Features")
    print(f"FPR = {num_retained_features} / {num_original_features}")
    print(f"Minimum FPR = {fpr:.4f}")

calculate_min_fpr()
<<<0.2857>>>