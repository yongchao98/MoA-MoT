def solve_fpr_problem():
    """
    Calculates the minimum Feature Retention Rate (FPR) based on the problem description.
    """
    print("Step 1: Define the original German semantic features for the concept 'mouldy'.")
    print("Let's assume the rich German concept needs to describe both cellar and fabric contexts.")
    print("- Core features (essential to any 'mouldy' smell): {'f_musty', 'f_damp'}")
    print("- Context-specific features: {'f_earthy' (for cellar), 'f_stale' (for fabric)}")
    print("The total set of original features is the union of these.")

    num_core_features = 2
    num_context_features = 2
    num_original_features = num_core_features + num_context_features

    print(f"\nTotal number of Original Features = {num_original_features}\n")

    print("Step 2: Determine the minimum features to retain in English.")
    print("The retained features must satisfy two conditions:")
    print(" A) Preserve the core meaning of 'mouldy'.")
    print(" B) Allow discrimination between 'mouldy cellar' and 'mouldy fabric'.")
    print("\nTo satisfy condition (A), we must retain at least 1 core feature (e.g., 'f_musty').")
    print("However, retaining only core features would map both 'mouldy cellar' and 'mouldy fabric' to the same set, failing condition (B).")
    print("To satisfy condition (B), we must also retain at least 1 context-specific feature (e.g., 'f_earthy').")
    print("\nTherefore, the minimum set includes 1 core feature + 1 context-specific feature.")

    num_retained_features = 2

    print(f"Minimum number of Retained Features = {num_retained_features}\n")

    print("Step 3: Calculate the minimum Feature Retention Rate (FPR).")
    print("FPR = Retained Features / Original Features")

    # Calculate the result
    fpr = num_retained_features / num_original_features

    # Print the final equation and result
    print(f"\nThe final calculation is:")
    print(f"{num_retained_features} / {num_original_features} = {fpr}")

solve_fpr_problem()
