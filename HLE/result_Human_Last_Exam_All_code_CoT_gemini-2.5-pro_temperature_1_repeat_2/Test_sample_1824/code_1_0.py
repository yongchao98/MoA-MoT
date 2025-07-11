import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """
    # Step 1: Model the German concepts with semantic features.
    # German has high granularity, so we define distinct features for each context.
    # Let's assume: f1='core_mouldy', f2='damp_earth', f3='stale_air', f4='decaying_textile'
    german_mouldy_cellar_features = {'core_mouldy', 'damp_earth', 'stale_air'}
    german_mouldy_fabric_features = {'core_mouldy', 'decaying_textile'}

    # Step 2: Determine the total set of original features.
    # This is the union of all features needed to describe both contexts in German.
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    # Step 3: Determine the minimum number of features to retain for discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric' in English, their
    # resulting feature sets must not be identical after mapping.
    
    # We must retain the 'core_mouldy' feature, otherwise the fundamental concept is lost.
    # To maintain distinction, we need to keep at least one other feature that is
    # present in one context but not the other.
    
    # Let's test retaining just two features: 'core_mouldy' and 'damp_earth'.
    retained_features_for_discrimination = {'core_mouldy', 'damp_earth'}
    num_retained_features = len(retained_features_for_discrimination)
    
    # Check if this minimal set works:
    # Map the German concepts to English by finding the intersection with the retained set.
    english_mouldy_cellar = german_mouldy_cellar_features.intersection(retained_features_for_discrimination)
    english_mouldy_fabric = german_mouldy_fabric_features.intersection(retained_features_for_discrimination)
    
    # Discrimination is successful if the sets are not equal.
    # english_mouldy_cellar -> {'core_mouldy', 'damp_earth'}
    # english_mouldy_fabric -> {'core_mouldy'}
    # Since they are different, discrimination is maintained with 2 features.
    # Retaining only 1 feature ('core_mouldy') would make both sets {'core_mouldy'}, failing the test.
    # Therefore, the minimum number of retained features is 2.

    # Step 4: Calculate the minimum FPR.
    min_fpr = num_retained_features / num_original_features

    # Print the explanation and the final equation.
    print("Problem Analysis:")
    print(f"1. Original German features for 'mouldy cellar': {german_mouldy_cellar_features}")
    print(f"2. Original German features for 'mouldy fabric': {german_mouldy_fabric_features}")
    print(f"3. Total unique original features (Union): {original_features}")
    print(f"   - Number of Original Features = {num_original_features}")
    print("\nFinding Minimum Retained Features for Discrimination:")
    print(f" - To distinguish the two concepts, their mapped feature sets in English must be different.")
    print(f" - We must retain 'core_mouldy'.")
    print(f" - We must also retain at least one distinguishing feature (e.g., 'damp_earth').")
    print(f" - Minimum features to retain: {retained_features_for_discrimination}")
    print(f"   - Minimum Number of Retained Features = {num_retained_features}")

    print("\nCalculation:")
    print(f"Minimum FPR = (Retained Features) / (Original Features)")
    print(f"Minimum FPR = {num_retained_features} / {num_original_features} = {min_fpr}")
    
    return min_fpr

# Execute the function and capture the final answer.
final_answer = calculate_min_fpr()
print(f"\n<<<__{final_answer}__>>>")
