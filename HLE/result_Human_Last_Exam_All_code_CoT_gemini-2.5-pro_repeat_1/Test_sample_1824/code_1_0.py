import sys

def calculate_minimum_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) for mapping a complex
    German olfactory concept to English while preserving discrimination.
    """

    # Step 1: Model the granular German concepts with a hypothetical set of features.
    # We assume 'mouldy cellar' and 'mouldy fabric' share core 'mouldy' features
    # but also have unique features reflecting German's lexical richness ('moderig', 'stickig', etc.).
    german_mouldy_cellar_features = {'fungal_decay', 'dampness', 'earthy', 'stagnant_air', 'cold_stone'}
    german_mouldy_fabric_features = {'fungal_decay', 'dampness', 'musty_textile', 'confined_space', 'organic_material'}

    # The total set of original features is the union of all unique features from both concepts.
    original_features = german_mouldy_cellar_features.union(german_mouldy_fabric_features)
    num_original_features = len(original_features)

    # Step 2: Determine the minimum set of features to retain in English.
    # To maintain the basic meaning of "mouldy", at least one core feature must be kept.
    # Let's assume 'fungal_decay' is the most essential core feature.
    core_retained = {'fungal_decay'}

    # To maintain discrimination between 'cellar' and 'fabric', we must retain
    # at least one unique feature for each concept that distinguishes them.
    # The most minimal way to do this is to pick one unique feature from each original set.
    # Let's pick 'earthy' for cellar and 'musty_textile' for fabric.
    distinguishing_retained = {'earthy', 'musty_textile'}

    # The minimum total set of retained features is the union of the core and distinguishing features.
    retained_features = core_retained.union(distinguishing_retained)
    num_retained_features = len(retained_features)

    # Step 3: Calculate the minimum FPR.
    if num_original_features == 0:
        # Avoid division by zero, though it's not expected in this scenario.
        min_fpr = 0
    else:
        min_fpr = num_retained_features / num_original_features

    # Step 4: Print the results clearly, showing the final equation.
    print(f"Total number of unique original German features: {num_original_features}")
    print(f"Original features considered: {original_features}\n")
    print(f"Minimum features to retain for discrimination: {num_retained_features}")
    print(f"Retained features considered: {retained_features}\n")
    print("The minimum Feature Retention Rate (FPR) is calculated as:")
    # The final print statement shows the full equation as requested.
    print(f"FPR = {num_retained_features} (Retained Features) / {num_original_features} (Original Features) = {min_fpr}")

    # The final answer for the platform.
    sys.stdout.write(f'<<<{min_fpr}>>>')

calculate_minimum_fpr()