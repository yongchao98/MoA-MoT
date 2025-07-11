import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """

    # Step 1: Model the granular German concepts with semantic features.
    # We define hypothetical feature sets. 'f_mould' is the core mouldy feature.
    # The other features represent the specific context (cellar vs. fabric).
    german_mouldy_cellar_features = {'f_mould', 'f_damp_earth', 'f_cold_air', 'f_stale'}
    german_mouldy_fabric_features = {'f_mould', 'f_damp_textile', 'f_chemical', 'f_dusty'}

    print("--- Step 1: Define German Concepts ---")
    print(f"Features for 'mouldy cellar': {german_mouldy_cellar_features}")
    print(f"Features for 'mouldy fabric': {german_mouldy_fabric_features}\n")

    # Step 2: Determine the total set of original features in German.
    # This is the union of all features from both concepts.
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    print("--- Step 2: Define Original Features (German) ---")
    print(f"Total set of original features is the union of the two sets.")
    print(f"Original Features: {original_features}")
    print(f"Total number of Original Features: {num_original_features}\n")

    # Step 3: Determine the minimum set of retained features in English.
    # To maintain discrimination, we need the core feature plus at least one
    # unique feature from each original concept.

    # Core features are the intersection of the two sets.
    core_features = german_mouldy_cellar_features.intersection(german_mouldy_fabric_features)

    # Unique features for each concept.
    unique_cellar_features = german_mouldy_cellar_features - german_mouldy_fabric_features
    unique_fabric_features = german_mouldy_fabric_features - german_mouldy_cellar_features
    
    # To minimize, we take one feature from each unique set.
    # We use list(set)[0] to pick an arbitrary but consistent feature.
    minimal_discriminator_cellar = {list(unique_cellar_features)[0]}
    minimal_discriminator_fabric = {list(unique_fabric_features)[0]}

    # The total set of retained features is the core feature(s) plus the minimal discriminators.
    retained_features = core_features.union(minimal_discriminator_cellar).union(minimal_discriminator_fabric)
    num_retained_features = len(retained_features)
    
    print("--- Step 3: Define Minimum Retained Features (English) ---")
    print(f"Core 'mouldy' feature(s) to retain: {core_features}")
    print(f"Minimum feature to distinguish 'cellar': {minimal_discriminator_cellar}")
    print(f"Minimum feature to distinguish 'fabric': {minimal_discriminator_fabric}")
    print(f"Total set of retained features in English map: {retained_features}")
    print(f"Total number of Retained Features: {num_retained_features}\n")
    
    # Step 4: Calculate the Feature Retention Rate (FPR).
    if num_original_features == 0:
        fpr = 0
    else:
        fpr = num_retained_features / num_original_features

    print("--- Step 4: Calculate Final FPR ---")
    print("FPR = Retained Features / Original Features")
    print(f"FPR = {num_retained_features} / {num_original_features}")
    print(f"The minimum FPR is: {fpr:.3f}")

    return fpr

# Run the calculation and store the result.
final_fpr = calculate_min_fpr()

# The final answer in the required format.
# The calculation is 3 / 7 = 0.42857...
# Rounding to three decimal places gives 0.429
print(f"\n<<<{final_fpr:.3f}>>>")
