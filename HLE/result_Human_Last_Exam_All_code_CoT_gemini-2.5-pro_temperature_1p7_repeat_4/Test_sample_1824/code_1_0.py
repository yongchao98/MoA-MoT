def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) based on a logical model
    of semantic features for olfactory concepts.
    """
    
    # Step 1 & 2: Define and calculate the minimal number of original features.
    # To have a shared concept ('mouldy'), there must be at least one common feature.
    min_common_features = 1
    # To have two distinct contexts ('cellar' vs 'fabric'), each must have at least one unique feature.
    min_cellar_specific_features = 1
    min_fabric_specific_features = 1

    # The total number of original features is the sum of these minimal components.
    total_original_features = min_common_features + min_cellar_specific_features + min_fabric_specific_features
    
    print(f"Modeling the original concepts in the high-granularity language (German):")
    print(f"- Minimum common features for 'mouldy': {min_common_features}")
    print(f"- Minimum specific features for 'cellar' context: {min_cellar_specific_features}")
    print(f"- Minimum specific features for 'fabric' context: {min_fabric_specific_features}")
    print(f"Total number of original features = {min_common_features} + {min_cellar_specific_features} + {min_fabric_specific_features} = {total_original_features}\n")

    # Step 3: Determine the minimum number of retained features.
    # To preserve the core 'mouldy' concept, we must retain at least one common feature.
    retained_common = 1
    # To maintain discrimination, we must also retain at least one specific feature.
    # Retaining just the common feature would make both concepts map to the same feature set.
    retained_specific = 1
    
    min_retained_features = retained_common + retained_specific
    
    print(f"Modeling the retained features for the low-granularity language (English):")
    print(f"- To preserve the core 'mouldy' concept, we need to retain {retained_common} common feature.")
    print(f"- To maintain discrimination, we need to retain {retained_specific} specific feature.")
    print(f"Minimum number of retained features = {retained_common} + {retained_specific} = {min_retained_features}\n")

    # Step 4: Calculate and print the minimum FPR.
    min_fpr_value = min_retained_features / total_original_features

    print("The final calculation for the minimum Feature Retention Rate (FPR) is:")
    print(f"FPR = Retained Features / Original Features")
    print(f"FPR = {min_retained_features} / {total_original_features}")
    print(f"The minimum FPR is: {min_fpr_value:.3f}")

calculate_min_fpr()