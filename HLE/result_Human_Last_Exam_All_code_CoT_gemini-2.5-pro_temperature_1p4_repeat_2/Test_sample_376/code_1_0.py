def calculate_clonal_expansion():
    """
    Calculates the clonal expansion score based on genetic and epigenetic events.
    """
    
    # Weights for the impact of each gene on clonal expansion.
    # Note: For tumor suppressors, a loss increases the score, so we use the positive magnitude of the weight.
    weight_oncogene_A = 0.5  # per additional copy
    weight_tsg_B = 0.7       # per lost copy
    weight_oncogene_C = 0.4  # per additional copy
    weight_tsg_D = 0.6       # per lost copy
    
    # --- Part 1: Calculate contributions from direct Copy Number Variations (CNVs) ---
    
    # Chromosome 1: gain of 3 copies, affecting Oncogene A.
    cnv_gain_chr1 = 3
    score_oncogene_A = cnv_gain_chr1 * weight_oncogene_A
    
    # Chromosome 2: loss of 2 copies, affecting Tumor Suppressor B.
    cnv_loss_chr2 = 2
    score_tsg_B = cnv_loss_chr2 * weight_tsg_B
    
    # Chromosome 3: gain of 2 copies, affecting Oncogene C.
    cnv_gain_chr3 = 2
    score_oncogene_C = cnv_gain_chr3 * weight_oncogene_C
    
    # --- Part 2: Calculate contributions from repressor overexpression ---
    
    # Repression is interpreted as a functional loss equivalent to losing 2 copies.
    functionally_lost_copies = 2
    
    # On Chromosome 1, repressor affects Tumor Suppressor D.
    score_tsg_D_repression = functionally_lost_copies * weight_tsg_D
    
    # On Chromosome 2, Tumor Suppressor B gene is already deleted,
    # so repression has no additional effect.
    
    # --- Part 3: Sum all contributions for the total score ---
    
    total_score = score_oncogene_A + score_tsg_B + score_oncogene_C + score_tsg_D_repression
    
    # Print the detailed equation and the final result.
    print("Clonal Expansion Score Calculation:")
    print(f"= (Oncogene A gain) + (Tumor Suppressor B loss) + (Oncogene C gain) + (Tumor Suppressor D repression)")
    print(f"= ({cnv_gain_chr1} * {weight_oncogene_A}) + ({cnv_loss_chr2} * {weight_tsg_B}) + ({cnv_gain_chr3} * {weight_oncogene_C}) + ({functionally_lost_copies} * {weight_tsg_D})")
    print(f"= {score_oncogene_A} + {score_tsg_B} + {score_oncogene_C} + {score_tsg_D_repression}")
    print(f"= {total_score}")

calculate_clonal_expansion()
<<<4.9>>>