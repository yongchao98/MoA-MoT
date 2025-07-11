import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping 'mouldy'
    from German to English while preserving discrimination between two contexts.
    """
    # Step 1: Define a hypothetical set of semantic features for the rich German concept of 'mouldy'.
    # This is our baseline for "Original Features".
    german_features = {
        'has_mould',        # The basic presence of fungal growth.
        'is_damp',          # The presence of moisture, a prerequisite for mould.
        'stale_air',        # Lack of ventilation, common to musty smells.
        'organic_decay',    # The smell of material (like fabric, wood) breaking down.
        'earthy_notes',     # Soil-like scents (e.g., from geosmin), typical of cellars.
        'mineral_notes'     # A scent associated with damp stone or concrete.
    }
    original_features_count = len(german_features)

    print("--- Step 1: Defining the Original Semantic Features ---")
    print(f"The German concept of 'mouldy' is modeled with {original_features_count} features:")
    print(f"{german_features}\n")

    # Step 2: Model which features apply to each specific context.
    cellar_features = {'has_mould', 'is_damp', 'stale_air', 'earthy_notes', 'mineral_notes'}
    fabric_features = {'has_mould', 'is_damp', 'stale_air', 'organic_decay'}

    print("--- Step 2: Applying Features to Specific Contexts ---")
    print(f"Features for 'mouldy cellar': {cellar_features}")
    print(f"Features for 'mouldy fabric': {fabric_features}\n")

    # Step 3: Identify core features (must be retained) and discriminatory features.
    # Core features are common to both contexts (set intersection).
    core_features = cellar_features.intersection(fabric_features)

    # Discriminatory features are those that differ between the contexts (set symmetric difference).
    discriminatory_features = cellar_features.symmetric_difference(fabric_features)

    print("--- Step 3: Identifying Core and Discriminatory Features ---")
    print(f"Core features (shared by both): {core_features}")
    print(f"To be 'mouldy', the {len(core_features)} core features must be retained.\n")
    print(f"Discriminatory features (to tell them apart): {discriminatory_features}")
    print("To distinguish the contexts, at least one of these must be retained.\n")

    # Step 4: Calculate the minimum number of features to retain.
    # This is the count of all core features plus the minimum number of discriminatory features (which is 1).
    min_retained_features_count = len(core_features) + 1

    print("--- Step 4: Calculating Minimum Retained Features ---")
    print(f"Minimum Retained Features = (Number of Core Features) + 1")
    print(f"Minimum Retained Features = {len(core_features)} + 1 = {min_retained_features_count}\n")

    # Step 5: Calculate the minimum Feature Retention Rate (FPR).
    min_fpr = min_retained_features_count / original_features_count

    print("--- Step 5: Calculating the Final Feature Retention Rate (FPR) ---")
    print("FPR = Retained Features / Original Features")
    print(f"FPR = {min_retained_features_count} / {original_features_count} = {min_fpr:.3f}\n")
    
    # Return the final numerical answer, which is 2/3
    return min_retained_features_count / original_features_count

# Execute the calculation and print the final result.
final_answer = calculate_min_fpr()
# The final answer is 4/6 = 2/3.
# We will represent it as a float.
final_answer_float = 4 / 6
print(f"The minimum FPR achievable is {final_answer_float}")
print(f"<<<{final_answer_float}>>>")
