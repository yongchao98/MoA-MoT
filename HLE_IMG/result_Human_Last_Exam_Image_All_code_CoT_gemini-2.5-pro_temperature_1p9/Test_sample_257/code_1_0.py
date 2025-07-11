def predict_nmr_features():
    """
    Analyzes the reaction of Pr-DAOTA with sulfuric acid and predicts the 1H NMR
    splitting pattern and integration for the most deshielded proton in the product.
    """

    # Step 1 & 2: Analyze the reaction and determine the product structure.
    # The reaction is a sulfonation. The N-propyl groups are strong activating ortho-, para-directors.
    # The most probable product, Compound 1, is formed by symmetrical di-sulfonation
    # at the positions ortho to the two nitrogen atoms. The original protons at these
    # positions are replaced by -SO3H groups.
    
    # Step 3: Identify the most deshielded proton in Compound 1.
    # After sulfonation, three types of aromatic protons remain:
    #   - Two equivalent protons ortho to the oxygen atom.
    #   - One unique proton on the central pyridine-like ring (let's call it H_c).
    # The proton H_c is part of a positively charged, pyridinium-like system and is flanked
    # by two nitrogen atoms. This highly electron-deficient environment makes it the most
    # deshielded proton, expected at a very high chemical shift (downfield).
    most_deshielded_proton = "the proton on the central ring (H_c)"

    # Step 4: Determine the splitting pattern and integration for this proton.
    
    # Splitting Pattern: A proton's signal is split by protons on adjacent carbon atoms.
    # The neighbors of the carbon holding H_c are both bridgehead carbons that do not
    # have any protons. With zero neighbors (n=0), the splitting pattern according to
    # the n+1 rule is a singlet.
    num_neighbors = 0
    splitting_pattern = "singlet"  # (0 + 1 = 1 peak)

    # Integration: The integration of a signal corresponds to the number of protons
    # that generate it. There is only one H_c proton in the molecule's structure.
    integration = 1

    # Print the final analysis
    print(f"The highest deshielded proton is {most_deshielded_proton}.")
    print("---")
    print("Prediction of its NMR signal:")
    print(f"Splitting Pattern: Based on {num_neighbors} neighboring protons, the signal is a {splitting_pattern}.")
    print(f"Integration: There is only one such proton, so the signal integrates to {integration}H.")
    print("\n--- Final Answer ---")
    print(f"The splitting pattern is a {splitting_pattern} and the integration is for {integration} proton.")


# Run the prediction
predict_nmr_features()