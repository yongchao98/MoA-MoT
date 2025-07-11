import math

def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) to discriminate
    between 'mouldy cellar' and 'mouldy fabric'.
    """

    # Step 1: Model the German concepts with semantic feature sets.
    # These include shared features ('mould', 'damp') and distinguishing features.
    german_mouldy_cellar_features = {'mould', 'damp', 'earthy', 'stale_air'}
    german_mouldy_fabric_features = {'mould', 'damp', 'organic_decay', 'sour'}

    # Step 2: Define "Original Features" as the union of all features involved.
    total_original_features_set = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(total_original_features_set)

    print("Step 1: Define Feature Sets")
    print(f"Features for 'mouldy cellar': {german_mouldy_cellar_features}")
    print(f"Features for 'mouldy fabric': {german_mouldy_fabric_features}")
    print("-" * 20)
    
    print("Step 2: Calculate Total Original Features")
    print(f"The total unique features across both concepts are: {total_original_features_set}")
    print(f"Total number of original features = {num_original_features}")
    print("-" * 20)

    # Step 3 & 4: Determine the minimum number of retained features for discrimination.
    # To maintain discrimination, the resulting feature sets after mapping must be different.
    # - If we retain a shared feature (e.g., 'damp'), both concepts map to {'damp'}, so they are not distinct.
    # - If we retain a single distinguishing feature (e.g., 'earthy'), 'cellar' maps to {'earthy'} and 'fabric' maps to {}.
    #   The sets are different, so discrimination is achieved.
    # Therefore, the minimum number of features we must retain is 1.
    min_retained_features = 1
    
    print("Step 3: Determine Minimum Retained Features for Discrimination")
    print("To distinguish the two concepts, the mapped feature sets cannot be identical.")
    print("Retaining just one unique feature (e.g., 'earthy') is sufficient.")
    print(f"Minimum number of retained features = {min_retained_features}")
    print("-" * 20)

    # Step 5: Calculate the minimum FPR.
    # FPR = Retained Features / Original Features
    min_fpr = min_retained_features / num_original_features

    print("Step 4: Calculate the Minimum Feature Retention Rate (FPR)")
    print("FPR = Minimum Retained Features / Total Original Features")
    # The final output prints each number in the equation as requested.
    print(f"Minimum FPR = {min_retained_features} / {num_original_features} = {min_fpr}")

if __name__ == "__main__":
    calculate_minimum_fpr()
