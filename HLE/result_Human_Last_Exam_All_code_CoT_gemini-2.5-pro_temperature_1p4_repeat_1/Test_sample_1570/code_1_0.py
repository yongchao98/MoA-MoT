import re

def calculate_mz_for_active_protein():
    """
    Calculates the m/z for a disulfide-linked peptide from XER22
    that would confirm the presence of the active therapeutic form.
    """
    protein_seq = "MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWRVVSSIEQKTEGAEKKQQMAREYREKIETELRDICNDVLSLLEKFLIPNACSQAESKVFYLKMKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWTSDRTQGCDEAEAGEGGEN"

    # Monoisotopic residue masses (mass of the amino acid minus H2O)
    aa_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    H2O_MASS = 18.01056
    H_ATOM_MASS = 1.007825
    PROTON_MASS = 1.007276

    # The second disulfide bridge is between Cys103 and Cys163.
    # We hypothesize a missed cleavage at K160.
    
    # Peptide containing Cys103
    peptide_C = "FLIPNACSQAESK"
    
    # Peptide containing Cys163, with one missed cleavage at K160
    peptide_D_missed = "LGLALNFSVFYYEILNSPEKACSLAK"
    
    def calculate_neutral_peptide_mass(sequence):
        """Calculates the neutral monoisotopic mass of a peptide sequence."""
        mass = H2O_MASS
        for aa in sequence:
            mass += aa_masses[aa]
        # round to 3rd decimal place as per instruction for intermediate steps
        return round(mass, 3)

    # Calculate masses of the individual peptides
    mass_pC = calculate_neutral_peptide_mass(peptide_C)
    mass_pD_missed = calculate_neutral_peptide_mass(peptide_D_missed)
    
    # Calculate mass of the disulfide-linked peptide pair
    # Mass = Mass(pC) + Mass(pD_missed) - 2 * Mass(H)
    linked_mass = mass_pC + mass_pD_missed - (2 * H_ATOM_MASS)
    linked_mass = round(linked_mass, 3)

    # Calculate m/z for different charge states
    print("Analysis for the second disulfide bridge with one missed cleavage:")
    print(f"Peptide 1: {peptide_C}")
    print(f"Peptide 2 (missed cleavage): {peptide_D_missed}")
    print("-" * 20)
    print(f"Mass of {peptide_C}: {mass_pC}")
    print(f"Mass of {peptide_D_missed}: {mass_pD_missed}")
    print(f"Mass of linked peptide pair: {mass_pC} + {mass_pD_missed} - 2 * {H_ATOM_MASS} = {linked_mass}")
    print("-" * 20)
    print("Potential m/z values:")

    for z in range(1, 6):
        mz = (linked_mass + z * PROTON_MASS) / z
        mz_rounded = round(mz, 3)
        print(f"  Charge (z) = +{z}: m/z = {mz_rounded}")

calculate_mz_for_active_protein()