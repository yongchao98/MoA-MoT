import sys

def solve_feature_retention():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving key distinctions.
    """
    # Step 1: Model the semantic features in German.
    # We hypothesize a set of features for the granular German concept of "mouldy"
    # that is rich enough to distinguish different contexts.
    # Let's assume the following 5 distinct features cover the space of 'mouldy cellar'
    # and 'mouldy fabric':
    # f1: Core 'mould' scent (essential for the general concept)
    # f2: 'Damp/earthy' scent (key differentiator for 'cellar')
    # f3: 'Stale/confined air' scent (nuance of 'cellar', can be lost)
    # f4: 'Decomposing organic matter' (common nuance in both, but can be lost if f1 is kept)
    # f5: 'Old/dusty textile' scent (key differentiator for 'fabric')
    # 'mouldy cellar' = {f1, f2, f3, f4}
    # 'mouldy fabric' = {f1, f4, f5}
    # The total pool of original features across these concepts is {f1, f2, f3, f4, f5}.
    original_features_count = 5

    # Step 2: Identify the minimum features to retain for discrimination in English.
    # To distinguish 'mouldy cellar' from 'mouldy fabric', we must retain:
    # 1. The core 'mould' concept (f1).
    # 2. The key feature for the 'cellar' context (f2: 'damp/earthy').
    # 3. The key feature for the 'fabric' context (f5: 'old textile').
    # The nuances (f3 and f4) can be lost while still meeting the minimum discrimination requirement.
    # Therefore, the minimum number of retained features is 3.
    retained_features_count = 3

    # Step 3: Calculate the Feature Retention Rate (FPR).
    # FPR = Retained Features / Original Features
    if original_features_count == 0:
        # Avoid division by zero, though not expected in this scenario.
        fpr = 0
    else:
        fpr = retained_features_count / original_features_count

    # Step 4: Print the final result including the equation.
    print(f"Hypothetical Original Features (German): {original_features_count}")
    print(f"Minimum Retained Features for Discrimination (English): {retained_features_count}")
    print("\nTo calculate the minimum Feature Retention Rate (FPR):")
    print(f"FPR = Retained Features / Original Features")
    print(f"FPR = {retained_features_count} / {original_features_count} = {fpr}")
    
    # Required for the final answer format
    sys.stdout.write(f"\n<<<{fpr}>>>")


solve_feature_retention()