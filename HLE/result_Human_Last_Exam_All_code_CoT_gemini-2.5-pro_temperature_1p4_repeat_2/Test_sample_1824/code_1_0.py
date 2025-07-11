import textwrap

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between contexts.
    """
    # Step 1: Model the German concepts with a plausible set of features.
    # We define a superset of features for the German concept 'moderig' (mouldy).
    original_features_german = {'f_musty', 'f_damp', 'f_earthy', 'f_decay', 'f_fungal'}
    num_original_features = len(original_features_german)

    print("Step-by-step analysis to find the Minimum Feature Retention Rate (FPR):")
    print("-" * 70)
    print("1. Defining the semantic feature sets based on the problem statement.")
    print(f"   - Original German 'moderig' features: {original_features_german}")
    print(f"   - Total number of original features: {num_original_features}\n")

    # Step 2: Model how different contexts activate subsets of these features in German.
    # 'mouldy cellar' might evoke dampness, earthiness, and fungal notes.
    activated_by_cellar = {'f_damp', 'f_earthy', 'f_fungal'}
    # 'mouldy fabric' might evoke mustiness, dampness, and decay (mildew).
    activated_by_fabric = {'f_musty', 'f_damp', 'f_decay'}

    print("2. Modeling contextual activation in German.")
    print(f"   - Features activated by 'cellar' context: {activated_by_cellar}")
    print(f"   - Features activated by 'fabric' context: {activated_by_fabric}\n")
    # In German, these sets are different, so 'moderiger keller' and 'moderiger stoff' are distinct.

    # Step 3: Identify core features that must be retained.
    # Core features are common to all contexts (intersection).
    core_features = activated_by_cellar.intersection(activated_by_fabric)
    
    explanation = textwrap.fill(
        "The problem states that 'core semantic features' must be preserved. We define these as the features common to both contexts.", 80
    )
    print("3. Identifying core features.")
    print(explanation)
    print(f"   - Core features to be retained = {core_features}\n")

    # Step 4: Determine the minimum features needed for discrimination in English.
    # The discrimination condition: (Retained ∩ Activated_Cellar) ≠ (Retained ∩ Activated_Fabric)
    
    # First, let's test if retaining only core features is enough.
    retained_test_set = core_features
    english_cellar_test = activated_by_cellar.intersection(retained_test_set)
    english_fabric_test = activated_by_fabric.intersection(retained_test_set)
    
    print("4. Checking if retaining only core features allows discrimination.")
    print(f"   - If we only retain {retained_test_set}:")
    print(f"     - English 'mouldy cellar' features: {english_cellar_test}")
    print(f"     - English 'mouldy fabric' features: {english_fabric_test}")
    print("   - The resulting sets are identical. Discrimination is lost.\n")

    # To regain discrimination, we must retain at least one feature that is not shared
    # between the two contexts. These are the features in the symmetric difference.
    differentiating_features = activated_by_cellar.symmetric_difference(activated_by_fabric)
    
    explanation2 = textwrap.fill(
        "To regain discrimination, we must add at least one feature from the set of features that differentiate the two contexts (their symmetric difference).", 80
    )
    print("5. Finding the minimum additional feature for discrimination.")
    print(explanation2)
    print(f"   - Differentiating features = {differentiating_features}")
    
    # The minimum set for retention is the core set plus one feature from the differentiating set.
    min_retained_features_set = core_features.union({list(differentiating_features)[0]})
    num_retained_features = len(min_retained_features_set)
    
    print(f"   - Therefore, the minimum retained set must contain the core feature(s) plus one differentiating feature.")
    print(f"   - Example minimum retained set: {min_retained_features_set}")
    print(f"   - Minimum number of retained features = {len(core_features)} (core) + 1 (differentiating) = {num_retained_features}\n")

    # Step 6: Calculate the minimum FPR.
    min_fpr = num_retained_features / num_original_features

    print("6. Calculating the final minimum FPR.")
    print("   Feature Retention Rate (FPR) = Retained Features / Original Features")
    print("   The equation is:")
    print(f"   FPR = {num_retained_features} / {num_original_features}")
    print(f"   Minimum FPR = {min_fpr}\n")
    
    return min_fpr

# Execute the function and print the final result.
final_fpr = calculate_min_fpr()
print("Final Answer:")
print(f"<<<{final_fpr}>>>")
