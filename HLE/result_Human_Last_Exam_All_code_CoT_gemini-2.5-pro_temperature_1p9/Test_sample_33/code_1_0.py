def analyze_protein_folding():
    """
    Analyzes FTIR data to explain the structural changes in a tardigrade protein
    during hydrogel formation.
    """
    
    # FTIR peak assignments for protein secondary structures
    peak_disordered = 1645  # cm^-1, broad, indicates random coil/disordered structures
    peak_alpha_helix = 1652   # cm^-1, sharp, indicates alpha-helical structures
    peak_beta_sheet_1 = 1618  # cm^-1, sharp, strong peak for anti-parallel beta-sheets
    peak_beta_sheet_2 = 1680  # cm^-1, sharp, weak peak for anti-parallel beta-sheets

    # Explanation based on experimental observations
    print("Step 1: Assigning FTIR peaks to secondary structures.")
    print(f"The initial protein is disordered, consistent with the broad peak at {peak_disordered} cm^-1.")
    print(f"The peak at {peak_alpha_helix} cm^-1 is characteristic of an alpha-helix.")
    print(f"The peaks at {peak_beta_sheet_1} cm^-1 and {peak_beta_sheet_2} cm^-1 are characteristic of beta-sheets.")
    print("\n" + "="*50 + "\n")

    print("Step 2: Interpreting the concentration titration experiment.")
    print("The experiment shows a dual increase in the peaks at:")
    print(f"- {peak_alpha_helix} cm^-1 (alpha-helix)")
    print(f"- {peak_beta_sheet_1} cm^-1 (beta-sheet)")
    print("This occurs as concentration increases, leading to gelation.")
    print("\n" + "="*50 + "\n")

    print("Conclusion:")
    print("The observation that signals for both alpha-helices and beta-sheets increase simultaneously")
    print("means that the initially disordered protein is folding into BOTH of these ordered structures")
    print("upon gelation. The correct explanation is that disordered structures fold into beta sheets and alpha helices.")

analyze_protein_folding()