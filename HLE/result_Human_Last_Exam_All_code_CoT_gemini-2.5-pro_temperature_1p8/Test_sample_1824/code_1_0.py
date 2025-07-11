def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while maintaining key distinctions.
    """
    # Step 1: Model the granular German concept of 'mouldy' ('moderig') with a set of features.
    # We hypothesize a set of 5 features for this nuanced concept.
    # Features: [f_musty, f_damp, f_fungal_growth, f_earthy, f_organic_decay]
    original_features_count = 5
    print(f"Assuming the original German concept has {original_features_count} semantic features.")
    print("Example Features: [f_musty, f_damp, f_fungal_growth, f_earthy, f_organic_decay]\n")

    # Step 2 & 3: Determine the minimum features to retain for discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric', we need:
    # 1. A core feature for 'mouldy' itself -> 'f_fungal_growth'
    # 2. A feature to identify the 'cellar' context -> 'f_earthy'
    # 3. A feature to identify the 'fabric' context -> 'f_organic_decay'
    # These three features are the absolute minimum required.
    retained_features_count = 3
    print(f"To distinguish 'mouldy cellar' from 'mouldy fabric', we must retain a minimum of {retained_features_count} features.")
    print("Retained Features: [f_fungal_growth, f_earthy, f_organic_decay]\n")

    # Step 4: Calculate the Feature Retention Rate (FPR).
    fpr = retained_features_count / original_features_count
    
    print("Calculating the minimum Feature Retention Rate (FPR)...")
    print("FPR = (Number of Retained Features) / (Number of Original Features)")
    print(f"FPR = {retained_features_count} / {original_features_count}")
    print(f"Minimum FPR = {fpr}")

calculate_minimum_fpr()
<<<0.6>>>