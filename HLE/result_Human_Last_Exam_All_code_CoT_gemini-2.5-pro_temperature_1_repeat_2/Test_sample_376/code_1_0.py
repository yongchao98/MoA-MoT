def calculate_clonal_expansion_score():
    """
    Calculates the clonal expansion score based on given CNVs and gene weights.
    """
    # Define copy number variations (CNVs)
    # Chromosome 1: gain of 3 copies
    gain_chr1 = 3
    # Chromosome 2: loss of 2 copies
    loss_chr2 = 2
    # Chromosome 3: gain of 2 copies
    gain_chr3 = 2

    # Define weights for gene impacts
    # Oncogenes (positive score for gains)
    weight_onc_A = 0.5
    weight_onc_C = 0.4
    # Tumor suppressors (positive score for losses, represented by negative weights)
    # We will use the absolute value for calculation
    weight_tsg_B = -0.7

    # --- Calculation ---
    # The score is only affected by oncogene gains and tumor suppressor losses.

    # 1. Chromosome 1: Gain of 3 copies affects Oncogene A.
    # Tumor suppressor D is also gained, but its weight is for loss, so its contribution is 0.
    score_onc_A = gain_chr1 * weight_onc_A

    # 2. Chromosome 2: Loss of 2 copies affects Tumor suppressor B.
    # The score increases with lost copies, so we use the absolute value of the weight.
    # Oncogene E is also lost, but its weight is for gain, so its contribution is 0.
    score_tsg_B = loss_chr2 * abs(weight_tsg_B)

    # 3. Chromosome 3: Gain of 2 copies affects Oncogene C.
    # Tumor suppressor F is also gained, but its weight is for loss, so its contribution is 0.
    score_onc_C = gain_chr3 * weight_onc_C

    # 4. Sum the scores for the final result.
    total_score = score_onc_A + score_tsg_B + score_onc_C

    # Print the detailed equation and the final score
    print("The clonal expansion score is calculated by summing the impacts of relevant gene copy number changes:")
    print("Score = (Oncogene A contribution) + (Tumor Suppressor B contribution) + (Oncogene C contribution)")
    print(f"Score = ({gain_chr1} copies * {weight_onc_A}) + ({loss_chr2} copies * {abs(weight_tsg_B)}) + ({gain_chr3} copies * {weight_onc_C})")
    print(f"Score = {score_onc_A} + {score_tsg_B} + {score_onc_C} = {total_score}")

calculate_clonal_expansion_score()
<<<3.7>>>