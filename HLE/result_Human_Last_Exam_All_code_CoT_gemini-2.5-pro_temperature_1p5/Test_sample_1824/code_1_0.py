import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping
    German 'mouldy' concepts to English while maintaining discrimination.
    """

    # Step 1: Model the German concepts with semantic features.
    # Let's define features:
    # f1: core feature for 'mouldy'
    # f2: nuance of 'earthy/damp soil' (typical for 'moderig', cellar)
    # f3: nuance of 'stale air/lack of ventilation' (typical for 'stickig', cellar)
    # f4: nuance of 'damp textile' (typical for 'muffig', fabric)

    german_mouldy_cellar_features = {'f1_core_mould', 'f2_earthy_soil', 'f3_stale_air'}
    german_mouldy_fabric_features = {'f1_core_mould', 'f4_damp_textile'}

    print("Step 1: Define features for the German concepts.")
    print(f"  - Features for 'mouldy cellar': {german_mouldy_cellar_features}")
    print(f"  - Features for 'mouldy fabric': {german_mouldy_fabric_features}")
    print("-" * 30)

    # Step 2: Determine the total set of original features for this task.
    # This is the union of all features involved in the concepts we need to discriminate.
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    print("Step 2: Determine the total number of original features.")
    print(f"  - The combined set of original features is: {original_features}")
    print(f"  - Total number of Original Features = {num_original_features}")
    print("-" * 30)


    # Step 3: Determine the minimum number of features that must be retained.
    # Condition A: The core semantic feature ('f1_core_mould') must be preserved.
    retained_features_count = 1  # for f1_core_mould

    # Condition B: To discriminate between 'cellar' and 'fabric', we need to retain
    # at least one of the features that makes them different.
    # The differentiating features are in the symmetric difference of the two sets.
    discriminating_features = german_mouldy_cellar_features.symmetric_difference(german_mouldy_fabric_features)

    # We must retain at least one of these discriminating features.
    retained_features_count += 1

    print("Step 3: Determine the minimum number of retained features.")
    print("  - To preserve the core meaning of 'mouldy', we must retain 'f1_core_mould'. (Retained Features = 1)")
    print(f"  - The features that differentiate the concepts are: {discriminating_features}")
    print("  - To maintain discrimination, we must retain at least one of these. (Retained Features = 1 + 1 = 2)")
    print(f"  - Minimum number of Retained Features = {retained_features_count}")
    print("-" * 30)

    # Step 4: Calculate the minimum FPR.
    min_fpr = retained_features_count / num_original_features

    print("Step 4: Calculate the Minimum Feature Retention Rate (FPR).")
    print("FPR = Retained Features / Original Features")
    # Final Answer Output
    print(f"Minimum FPR = {retained_features_count} / {num_original_features} = {min_fpr}")

if __name__ == "__main__":
    calculate_min_fpr()