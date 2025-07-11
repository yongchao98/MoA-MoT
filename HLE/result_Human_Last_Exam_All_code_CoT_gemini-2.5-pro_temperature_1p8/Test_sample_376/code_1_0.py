def calculate_clonal_expansion_score():
    """
    Calculates the clonal expansion score for a tumor based on given genetic data.
    """

    # --- Data Definition ---
    # Copy number variations (CNVs)
    cnv_chr1 = 3  # Gain of 3 copies
    cnv_chr2 = -2 # Loss of 2 copies
    cnv_chr3 = 2  # Gain of 2 copies

    # Gene impact weights
    weights = {
        'Oncogene A': 0.5,
        'Tumor suppressor B': -0.7,
        'Oncogene C': 0.4,
        'Tumor suppressor D': -0.6,
        'Oncogene E': 0.3,
        'Tumor suppressor F': -0.5,
    }

    # Repressor status
    repressor_on_chr1 = True
    repressor_on_chr2 = True
    repressor_on_chr3 = False
    
    print("Calculating the clonal expansion score:\n")

    # --- Calculations ---
    # Chromosome 1
    impact_oncogene_A = cnv_chr1 * weights['Oncogene A']
    impact_ts_D_raw = cnv_chr1 * weights['Tumor suppressor D']
    # The repressor negates the benefit of the tumor suppressor gain
    impact_ts_D = 0 if repressor_on_chr1 else impact_ts_D_raw
    
    # Chromosome 2
    impact_ts_B = cnv_chr2 * weights['Tumor suppressor B']
    impact_oncogene_E = cnv_chr2 * weights['Oncogene E']

    # Chromosome 3
    impact_oncogene_C = cnv_chr3 * weights['Oncogene C']
    impact_ts_F = cnv_chr3 * weights['Tumor suppressor F']
    
    # Store all final impacts for the final equation
    all_impacts = [
        impact_oncogene_A,
        impact_ts_D,
        impact_ts_B,
        impact_oncogene_E,
        impact_oncogene_C,
        impact_ts_F
    ]
    
    # --- Output ---
    print(f"Oncogene A (Chr 1 Gain): {cnv_chr1} copies * {weights['Oncogene A']}/copy = {impact_oncogene_A:.1f}")
    print(f"Tumor Suppressor D (Chr 1 Gain): Impact is {impact_ts_D_raw:.1f} but negated by repressor -> {impact_ts_D:.1f}")
    print(f"Tumor Suppressor B (Chr 2 Loss): {cnv_chr2} copies * {weights['Tumor suppressor B']}/copy = {impact_ts_B:.1f}")
    print(f"Oncogene E (Chr 2 Loss): {cnv_chr2} copies * {weights['Oncogene E']}/copy = {impact_oncogene_E:.1f}")
    print(f"Oncogene C (Chr 3 Gain): {cnv_chr3} copies * {weights['Oncogene C']}/copy = {impact_oncogene_C:.1f}")
    print(f"Tumor Suppressor F (Chr 3 Gain): {cnv_chr3} copies * {weights['Tumor suppressor F']}/copy = {impact_ts_F:.1f}")
    
    # Constructing the final equation string with all numbers
    equation_str = " + ".join([f"({val:.1f})" if val < 0 else f"{val:.1f}" for val in all_impacts])
    
    # Calculate the total score
    total_score = sum(all_impacts)

    print("\n---")
    print("Final Equation:")
    print(f"{equation_str}")
    
    print("\nTotal Clonal Expansion Score:")
    print(f"{total_score:.1f}")

calculate_clonal_expansion_score()
<<<2.1>>>