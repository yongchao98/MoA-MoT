def peptide_synthesis_analysis():
    """
    Analyzes the synthesis of a 100aa peptide with an unnatural amino acid
    and calculates its theoretical mass.
    """
    
    # --- Part 1: Explanation of the Synthesis Strategy ---

    print("--- Recommended Synthesis Strategy ---")
    print("For a 100aa peptide containing an unnatural amino acid like azido phenylalanine, the most helpful technique is Native Chemical Ligation (NCL).")
    print("\nThis strategy involves:")
    print("1. Synthesizing two smaller, manageable peptide fragments (~50 aa each) using Solid-Phase Peptide Synthesis (SPPS).")
    print("2. Incorporating the unnatural amino acid (X) into one of the fragments during its synthesis.")
    print("3. Chemically ligating the two fragments at a Cysteine residue to form the full-length peptide.\n")

    # --- Part 2: Hypothetical Peptide and Mass Calculation ---

    # Define a hypothetical 100-amino-acid peptide for demonstration.
    # It has a Cysteine (C) at position 50, ideal for NCL,
    # and Azido Phenylalanine (X) at position 52.
    # Sequence: M + A(48) + C + A + X + A(48)
    peptide_sequence = "M" + "A" * 48 + "C" + "A" + "X" + "A" * 48
    
    cys_position = peptide_sequence.find('C')
    fragment1_seq = peptide_sequence[:cys_position]
    fragment2_seq = peptide_sequence[cys_position:]

    print("--- Illustrative Example for NCL ---")
    print(f"Hypothetical 100aa Sequence: {peptide_sequence[:25]}...{peptide_sequence[75:]}")
    print(f"Fragment 1 ({len(fragment1_seq)} aa) to synthesize with C-terminal thioester: {fragment1_seq[:25]}...")
    print(f"Fragment 2 ({len(fragment2_seq)} aa) to synthesize with N-terminal Cysteine: {fragment2_seq[:25]}...\n")

    # Monoisotopic masses of amino acid residues (mass of free aa - mass of H2O)
    residue_masses = {
        'A': 71.03711,   # Alanine
        'C': 103.00919,  # Cysteine
        'M': 131.04049,  # Methionine
        # X = p-azido-L-phenylalanine (AzF) residue. Formula: C9H8N4O
        'X': 188.06981,
    }
    # Mass of a water molecule (H2O)
    water_mass = 18.010565

    print("--- Theoretical Mass Calculation ---")
    total_residue_mass = 0
    mass_components = []

    # Calculate the sum of all residue masses
    for aa in peptide_sequence:
        mass = residue_masses[aa]
        total_residue_mass += mass
        mass_components.append(str(mass))

    # Final peptide mass = (Sum of residue masses) + (Mass of one water molecule for terminals)
    final_mass = total_residue_mass + water_mass
    
    print("The final mass is the sum of all constituent residue masses plus one water molecule (for the terminal -H and -OH groups).")
    print("Equation (showing first and last 3 residue masses):")
    equation_str = " + ".join(mass_components[:3]) + " + ... + " + " + ".join(mass_components[-3:])
    print(f"Sum of Residue Masses = {equation_str}")
    print(f"Total Sum of Residue Masses = {total_residue_mass:.5f} Da")
    print(f"Final Peptide Mass = (Total Sum of Residue Masses) + Mass(H2O)")
    print(f"Final Peptide Mass = {total_residue_mass:.5f} + {water_mass:.5f} = {final_mass:.5f} Da")

# Execute the analysis and calculation
peptide_synthesis_analysis()
