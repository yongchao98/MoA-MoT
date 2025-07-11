def calculate_peptide_mass():
    """
    Calculates and prints the monoisotopic molecular weight of a peptide.

    The script defines a representative 100-amino-acid peptide sequence
    containing azido phenylalanine (X). It then calculates the total mass
    by summing the masses of the individual amino acid residues and adding
    the mass of a water molecule for the termini.
    """
    # Monoisotopic masses of amino acid residues (mass of AA - mass of H2O)
    # 'X' is p-azido-L-phenylalanine (C9H8N4O residue)
    residue_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
        'X': 188.06981  # p-azido-L-phenylalanine residue mass
    }
    mass_of_water = 18.01528 # To add for the terminal H and OH groups

    # A representative 100aa sequence based on the user's query M...KAVCLXVIGATR...A
    # Sequence: M + 39*A + KAVCLXVIGATR + 47*A + A = 100 aa
    peptide_sequence = "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKAVCLXVIGATRAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    print(f"Calculating molecular weight for sequence: {peptide_sequence}\n")

    total_mass = 0.0
    mass_components = []

    # Sum the masses of all residues in the sequence
    for aa in peptide_sequence:
        if aa in residue_masses:
            mass = residue_masses[aa]
            total_mass += mass
            mass_components.append(str(mass))
        else:
            print(f"Warning: Amino acid '{aa}' not found in mass dictionary. Skipping.")
            return

    # Add the mass of a water molecule for the N- and C-termini
    total_mass += mass_of_water
    mass_components.append(str(mass_of_water))

    # Create the full equation string
    equation_str = " + ".join(mass_components)

    print("Final Equation:")
    # The equation is very long, so we print it in a more readable format if needed,
    # but here we will print it as one line to meet the prompt's requirement.
    print(f"{equation_str} = {total_mass:.5f}")

    print("\n--------------------------------------------------")
    print(f"Peptide Sequence Length: {len(peptide_sequence)} aa")
    print(f"Monoisotopic Molecular Weight: {total_mass:.5f} Da")
    print("--------------------------------------------------")

if __name__ == '__main__':
    calculate_peptide_mass()