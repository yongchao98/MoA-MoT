import math

def calculate_min_fpr():
    """
    Calculates the minimum Feature Retention Rate (FPR) based on a logical model
    of semantic features for German olfactory concepts.
    """
    # Step 1: Model the original German concepts.
    # We assume the most minimal logical structure for a concept that has both
    # a "core" meaning and a "nuance", as described in the problem.
    # - 1 core feature (e.g., 'musty'), common to both concepts.
    # - 1 nuance feature (e.g., 'cellar-specific'), unique to each concept.
    num_core_features = 1
    num_nuance_features = 1

    # The total number of features in the original German concept.
    original_features = num_core_features + num_nuance_features

    # Step 2: Determine the minimum features to retain for discrimination.
    # To distinguish 'mouldy cellar' from 'mouldy fabric' in English,
    # their retained feature sets must be different. The most minimal way to ensure this
    # is to retain only the single, unique nuance feature for each concept.
    # Retaining the core feature for both would make them identical.
    retained_features = 1

    # Step 3: Calculate the minimum Feature Retention Rate (FPR).
    # FPR = Retained Features / Original Features
    min_fpr = retained_features / original_features

    # Print the explanation and the final equation.
    print("### Step-by-Step Calculation ###")
    print(f"1. Model Original Concept: Each German concept ('mouldy cellar', 'mouldy fabric') is modeled with:")
    print(f"   - {num_core_features} core feature (e.g., 'musty')")
    print(f"   - {num_nuance_features} nuance feature (e.g., 'earthy')")
    print(f"   Total Original Features per concept = {original_features}")
    print("\n2. Determine Minimum Retained Features for Discrimination:")
    print(f"   To tell the concepts apart, we only need to keep the single distinguishing nuance feature.")
    print(f"   Minimum Retained Features per concept = {retained_features}")
    print("\n3. Calculate Minimum Feature Retention Rate (FPR):")
    print(f"   FPR = Retained Features / Original Features")
    # The final equation shows each number used in the calculation.
    print("\nFinal Equation:")
    print(f"{retained_features} / {original_features} = {min_fpr}")

calculate_min_fpr()

# The final numerical answer is 0.5
print("\n<<<0.5>>>")