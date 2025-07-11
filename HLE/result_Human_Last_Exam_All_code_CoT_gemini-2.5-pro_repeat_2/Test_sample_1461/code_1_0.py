def predict_helix_pattern():
    """
    Predicts the helical pattern for an alternating alpha/epsilon-peptide foldamer
    based on the established trend for alpha/X-peptide helices.
    """

    # Database of backbone lengths for X-amino acids
    backbone_atoms = {
        'beta': 4,
        'gamma': 5,
        'delta': 6,
        'epsilon': 7
    }

    # Base case from known experimental data: alpha/beta-peptides
    base_monomer_type = 'beta'
    base_atoms = backbone_atoms[base_monomer_type]
    base_helix_small_ring = 10
    base_helix_large_ring = 12

    # Target monomer from the problem
    target_monomer_type = 'epsilon'
    target_atoms = backbone_atoms[target_monomer_type]

    # The trend shows that for each additional atom in the X-monomer backbone,
    # the helical ring sizes both increase by 1.
    atom_difference = target_atoms - base_atoms

    # Extrapolate to find the new helical pattern
    predicted_small_ring = base_helix_small_ring + atom_difference
    predicted_large_ring = base_helix_large_ring + atom_difference

    print("--- Prediction of Helical Pattern for an alpha/epsilon-Peptide ---")
    print(f"Base case: alpha/{base_monomer_type}-peptides ({base_atoms} backbone atoms) form a {base_helix_small_ring}/{base_helix_large_ring}-helix.")
    print(f"Target monomer: epsilon-amino acid has {target_atoms} backbone atoms.")
    print(f"Difference in backbone atoms compared to base case: {target_atoms} - {base_atoms} = {atom_difference}")
    print("\nApplying the established linear trend:")
    print(f"Predicted small ring size: {base_helix_small_ring} + {atom_difference} = {predicted_small_ring}")
    print(f"Predicted large ring size: {base_helix_large_ring} + {atom_difference} = {predicted_large_ring}")
    print("\n-----------------------------------------------------------------")
    print(f"The most likely helical pattern is the {predicted_small_ring}/{predicted_large_ring}-helix.")
    print("-----------------------------------------------------------------")

predict_helix_pattern()
<<<B>>>