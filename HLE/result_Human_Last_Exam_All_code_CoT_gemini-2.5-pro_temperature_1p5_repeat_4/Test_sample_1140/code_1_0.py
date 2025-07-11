def calculate_peptide_fragment_mw():
    """
    Calculates the molecular weight of a peptide fragment and prints the detailed equation.
    """
    # Average molecular weights of amino acids (in g/mol)
    amino_acid_mw = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10,
        'C': 121.16, 'Q': 146.14, 'E': 147.13, 'G': 75.07,
        'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
        'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
        'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
        'X': 206.21  # Special case for p-Azido-L-phenylalanine
    }
    water_mw = 18.015

    sequence = "KAVCLXVIGATR"
    
    print(f"Calculating molecular weight for the peptide fragment: {sequence}\n")
    
    # --- Step 1: Sum of individual amino acid masses ---
    total_aa_mass = 0
    mass_equation_parts = []
    
    for aa in sequence:
        mass = amino_acid_mw.get(aa.upper())
        if mass is None:
            print(f"Error: Amino acid '{aa}' not found in the library.")
            return
        total_aa_mass += mass
        mass_equation_parts.append(str(mass))
        
    print("Step 1: Sum of individual amino acid masses")
    sum_equation = " + ".join(mass_equation_parts)
    print(f"{sum_equation} = {total_aa_mass:.2f} g/mol")
    print("-" * 30)

    # --- Step 2: Calculate mass of water lost ---
    num_peptide_bonds = len(sequence) - 1
    total_water_mass = num_peptide_bonds * water_mw
    
    print("Step 2: Calculate mass of water lost from peptide bond formation")
    print(f"Number of peptide bonds: {num_peptide_bonds}")
    print(f"Mass of water (H2O): {water_mw} g/mol")
    print(f"{num_peptide_bonds} * {water_mw} = {total_water_mass:.3f} g/mol")
    print("-" * 30)

    # --- Step 3: Final molecular weight ---
    final_mw = total_aa_mass - total_water_mass
    
    print("Step 3: Calculate final molecular weight of the fragment")
    print("(Total AA Mass) - (Total Water Mass)")
    print(f"{total_aa_mass:.2f} - {total_water_mass:.3f} = {final_mw:.3f} g/mol\n")
    print(f"Final Molecular Weight of '{sequence}': {final_mw:.3f} g/mol")

# Execute the calculation
calculate_peptide_fragment_mw()