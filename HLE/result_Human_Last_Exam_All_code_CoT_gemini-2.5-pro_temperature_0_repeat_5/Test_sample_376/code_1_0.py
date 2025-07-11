def calculate_clonal_expansion():
    """
    Calculates the clonal expansion score based on given genetic data and weights.
    """
    # --- Data from the problem statement ---

    # Copy number changes for each gene
    # Chromosome 1: Gain of 3 copies
    oncogene_A_change = 3
    tumor_suppressor_D_change = 3
    # Chromosome 2: Loss of 2 copies
    tumor_suppressor_B_change = -2
    oncogene_E_change = -2
    # Chromosome 3: Gain of 2 copies
    tumor_suppressor_F_change = 2
    oncogene_C_change = 2

    # Weights for impact. We use the absolute value as loss of a tumor suppressor
    # contributes positively to clonal expansion.
    weights = {
        'Oncogene A': 0.5,       # per additional copy
        'Tumor suppressor B': 0.7, # per lost copy
        'Oncogene C': 0.4,       # per additional copy
        'Tumor suppressor D': -0.6, # per lost copy (no weight for gain)
        'Oncogene E': 0.3,       # per additional copy (no weight for loss)
        'Tumor suppressor F': -0.5  # per lost copy (no weight for gain)
    }

    # --- Calculation ---

    # Calculate contribution from Oncogene A (gain)
    score_A = 0
    if oncogene_A_change > 0:
        score_A = oncogene_A_change * weights['Oncogene A']

    # Calculate contribution from Tumor suppressor B (loss)
    score_B = 0
    if tumor_suppressor_B_change < 0:
        # Number of lost copies is the absolute value of the change
        num_lost = abs(tumor_suppressor_B_change)
        score_B = num_lost * weights['Tumor suppressor B']

    # Calculate contribution from Oncogene C (gain)
    score_C = 0
    if oncogene_C_change > 0:
        score_C = oncogene_C_change * weights['Oncogene C']

    # Other genes do not contribute as the type of change does not have a specified weight
    # Tumor suppressor D: Gained, but weight is for loss. Contribution = 0.
    # Oncogene E: Lost, but weight is for gain. Contribution = 0.
    # Tumor suppressor F: Gained, but weight is for loss. Contribution = 0.
    score_D = 0
    score_E = 0
    score_F = 0
    
    # --- Output ---

    total_score = score_A + score_B + score_C + score_D + score_E + score_F

    print("Clonal Expansion Score Calculation:")
    print(f"1. Oncogene A (gain of {oncogene_A_change}): {oncogene_A_change} copies * {weights['Oncogene A']}/copy = {score_A}")
    print(f"2. Tumor Suppressor B (loss of {abs(tumor_suppressor_B_change)}): {abs(tumor_suppressor_B_change)} copies * {weights['Tumor suppressor B']}/copy = {score_B}")
    print(f"3. Oncogene C (gain of {oncogene_C_change}): {oncogene_C_change} copies * {weights['Oncogene C']}/copy = {score_C}")
    print("4. Other genes (D, E, F) have changes for which no weights are provided, so their contribution is 0.")
    
    print("\nFinal Equation:")
    # We show the full equation including the zero-contribution terms for clarity
    print(f"{score_A} (from A) + {score_B} (from B) + {score_C} (from C) + {score_D} (from D) + {score_E} (from E) + {score_F} (from F) = {total_score}")
    
    print(f"\nThe final clonal expansion score is: {total_score}")

calculate_clonal_expansion()