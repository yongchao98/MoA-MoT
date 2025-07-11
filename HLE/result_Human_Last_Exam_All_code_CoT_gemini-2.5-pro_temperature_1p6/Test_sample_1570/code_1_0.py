import re

def solve_protein_mz():
    """
    This function calculates the m/z values for disulfide-linked tryptic peptides
    of protein XER22 to identify evidence of its active form.
    """
    
    # Monoisotopic residue masses and other constants
    residue_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00918, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04048, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    MASS_H2O = 18.01056
    MASS_H = 1.007825
    MASS_PROTON = 1.007276

    full_sequence = "MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWRVVSSIEQKTEGAEKKQQMAREYREKIETELRDICNDVLSLLEKFLIPNACSQAESKVFYLKMKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWTSDRTQGCDEAEAGEGGEN"

    # Step 1 & 2: Identify the tryptic peptides containing the cysteines for the active bridges.
    # This is done based on the analysis of trypsin cleavage sites (after K or R, not before P).
    # From manual analysis:
    # Bridge 1 involves linking Cys30 and Cys220.
    peptide_1a = "LAEQAERYDDMAACMK" # Contains Cys from "MAACM" region
    peptide_1b = "TQGCDEAEAGEGGEN"  # Contains Cys from "TQGCDEAEAGEG" region

    # Bridge 2 involves linking Cys133 and Cys175.
    peptide_2a = "FLIPNACSQAESK"   # Contains Cys from "NACSQAESK" region
    peptide_2b = "ACSLAK"          # Contains Cys from "PEKACSLAKTAFDEA" region

    peptides = {
        "Bridge 1": (peptide_1a, peptide_1b),
        "Bridge 2": (peptide_2a, peptide_2b)
    }
    
    print("--- Calculating m/z values for disulfide-linked peptides ---")

    for bridge_name, (p1_seq, p2_seq) in peptides.items():
        # Step 3: Calculate mass of the linked pair
        
        # Calculate neutral mass of peptide 1
        p1_res_mass = sum(residue_masses[aa] for aa in p1_seq)
        p1_neutral_mass = p1_res_mass + MASS_H2O
        
        # Calculate neutral mass of peptide 2
        p2_res_mass = sum(residue_masses[aa] for aa in p2_seq)
        p2_neutral_mass = p2_res_mass + MASS_H2O

        # Calculate mass of the disulfide-linked complex
        linked_mass = p1_neutral_mass + p2_neutral_mass - (2 * MASS_H)

        print(f"\nAnalysis for {bridge_name}:")
        print(f"Peptide 1: {p1_seq}")
        print(f"Peptide 2: {p2_seq}")
        print(f"Combined Neutral Mass [M]: {linked_mass:.3f} Da")

        # Step 4: Calculate m/z values for different charge states
        print("Potential m/z values [M+zH]^(z+):")
        for z in range(2, 5): # Common charge states for these sizes
            mz_value = (linked_mass + (z * MASS_PROTON)) / z
            print(f"  z={z}: {mz_value:.3f}")

    print("\n--- Comparing with options ---")
    print("Bridge 1 calculated m/z values (e.g., 1103.440 for z=3) do not match any options.")
    print("Bridge 2 calculated m/z values (e.g., 998.995 for z=2) do not match any options.")
    print("Therefore, assuming no missed cleavages or modifications, none of the direct calculations match the provided choices.")

solve_protein_mz()