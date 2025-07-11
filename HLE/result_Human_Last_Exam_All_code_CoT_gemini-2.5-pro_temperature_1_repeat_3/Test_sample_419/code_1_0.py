def calculate_peptide_mw():
    """
    Calculates the molecular weight of the peptide GSTA.
    The MUC1 protein contains many repeats of a similar sequence.
    """
    # Average molecular weights of amino acids (in g/mol)
    amino_acid_weights = {
        'G': 75.07,   # Glycine
        'S': 105.09,  # Serine
        'T': 119.12,  # Threonine
        'A': 89.09    # Alanine
    }

    # Weight of a water molecule (H2O) to be subtracted for each peptide bond
    water_weight = 18.015

    peptide_sequence = "GSTA"
    num_amino_acids = len(peptide_sequence)
    num_peptide_bonds = num_amino_acids - 1

    # Sum the weights of the individual amino acids
    sum_of_residue_weights = sum(amino_acid_weights[aa] for aa in peptide_sequence)
    
    # Calculate the total weight of water molecules removed
    total_water_weight_removed = num_peptide_bonds * water_weight
    
    # Final peptide molecular weight
    peptide_mw = sum_of_residue_weights - total_water_weight_removed
    
    print(f"Calculating molecular weight for peptide: {peptide_sequence}")
    print("-" * 40)
    print(f"Sum of residue weights (G+S+T+A):")
    g_weight = amino_acid_weights['G']
    s_weight = amino_acid_weights['S']
    t_weight = amino_acid_weights['T']
    a_weight = amino_acid_weights['A']
    print(f"{g_weight} + {s_weight} + {t_weight} + {a_weight} = {sum_of_residue_weights:.2f} g/mol")
    
    print(f"\nNumber of peptide bonds: {num_peptide_bonds}")
    print(f"Weight of water removed per bond: {water_weight} g/mol")
    print(f"Total weight removed: {num_peptide_bonds} * {water_weight} = {total_water_weight_removed:.3f} g/mol")

    print("\nFinal Equation:")
    print(f"({g_weight} + {s_weight} + {t_weight} + {a_weight}) - ({num_peptide_bonds} * {water_weight}) = {peptide_mw:.3f} g/mol")
    
    print("-" * 40)
    print(f"The calculated molecular weight of GSTA is: {peptide_mw:.3f} g/mol")

calculate_peptide_mw()