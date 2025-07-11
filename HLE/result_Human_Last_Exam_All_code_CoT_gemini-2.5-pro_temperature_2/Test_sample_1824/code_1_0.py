def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping
    German 'mouldy' concepts to English while maintaining discrimination.
    """

    # Step 1: Define feature sets for the German concepts based on the problem description.
    # We need shared features for 'mouldy' and unique features for 'cellar' and 'fabric'.
    german_mouldy_cellar_features = {'f_mould', 'f_damp', 'f_earthy'}
    german_mouldy_fabric_features = {'f_mould', 'f_damp', 'f_organic_decay'}

    # Step 2: The total set of "Original Features" is the union of all features
    # required to describe both concepts in German.
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    # Step 3: Determine the minimum number of features that must be retained in English.
    # Condition A: The core concept of 'mouldy' must be preserved. We assume this means 'f_mould' must be retained.
    # Condition B: Basic discrimination must be possible. To make the two concepts different,
    # the mapping only needs to retain at least one of the features that makes them unique ('f_earthy' or 'f_organic_decay').
    # Therefore, the smallest possible set of retained features would be {'f_mould', 'f_earthy'} or {'f_mould', 'f_organic_decay'}.
    # The size of this minimum set is 2.
    num_retained_features = 2

    # Step 4: Calculate the minimum FPR.
    min_fpr = num_retained_features / num_original_features

    # Print the explanation and the calculation
    print("Plan:")
    print("1. Define semantic features for 'mouldy cellar' and 'mouldy fabric'.")
    print("   - 'mouldy cellar': " + str(german_mouldy_cellar_features))
    print("   - 'mouldy fabric': " + str(german_mouldy_fabric_features))
    print("\n2. The total set of original features is the union of both sets:")
    print("   - Original Features: " + str(original_features))
    print(f"   - Number of Original Features = {num_original_features}")
    print("\n3. Determine the minimum features to retain in the English mapping.")
    print("   - To preserve the core meaning, we must keep 'f_mould'.")
    print("   - To tell 'cellar' and 'fabric' apart, we only need to keep one of the distinguishing features (e.g., 'f_earthy').")
    print("   - A minimal set is {'f_mould', 'f_earthy'}.")
    print(f"   - Minimum Number of Retained Features = {num_retained_features}")
    print("\n4. Calculate the minimum Feature Retention Rate (FPR).")
    print(f"   FPR = Retained Features / Original Features")
    print(f"\nFinal Calculation: {num_retained_features} / {num_original_features} = {min_fpr}")
    
    # Final answer in the required format
    print(f"\n<<<{min_fpr}>>>")

calculate_min_fpr()