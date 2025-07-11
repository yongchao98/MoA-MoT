def calculate_clonal_expansion_score():
    """
    Calculates the clonal expansion score based on Copy Number Variations (CNVs)
    and gene-specific weights.
    """
    # Chromosome 1: Gain of 3 copies affecting Oncogene A
    gain_oncogene_A = 3
    weight_oncogene_A = 0.5
    score_A = gain_oncogene_A * weight_oncogene_A

    # Chromosome 2: Loss of 2 copies affecting Tumor Suppressor B
    # The loss of a tumor suppressor increases the clonal expansion score.
    # We multiply by -1 to make the final contribution positive.
    loss_tsg_B = 2
    weight_tsg_B = -0.7
    score_B = -1 * loss_tsg_B * weight_tsg_B

    # Chromosome 3: Gain of 2 copies affecting Oncogene C
    gain_oncogene_C = 2
    weight_oncogene_C = 0.4
    score_C = gain_oncogene_C * weight_oncogene_C

    # Total score is the sum of the individual contributions
    total_score = score_A + score_B + score_C
    
    # Print the equation with all the numbers used in the calculation
    # The calculation for the tumor suppressor B is shown as (-2 * -0.7) to clearly represent the loss effect
    print(f"Clonal Expansion Score Calculation:")
    print(f"({gain_oncogene_A} * {weight_oncogene_A}) + (-{loss_tsg_B} * {weight_tsg_B}) + ({gain_oncogene_C} * {weight_oncogene_C}) = {total_score}")

calculate_clonal_expansion_score()
<<<3.7>>>