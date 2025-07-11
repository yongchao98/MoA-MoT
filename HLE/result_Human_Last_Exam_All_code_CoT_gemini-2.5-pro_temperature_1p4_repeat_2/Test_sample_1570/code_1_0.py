import sys

def calculate_mz():
    """
    Identifies tryptic peptides involved in disulfide bridges of protein XER22,
    calculates their combined mass, and determines the m/z value that matches
    the experimental options, assuming a sequence typo and a sodium adduct.
    """
    # Monoisotopic residue masses of amino acids (mass of AA - mass of H2O)
    aa_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    # Mass of constituents
    H2O_MASS = 18.010565
    H_MASS = 1.007825
    PROTON_MASS = 1.007276
    NA_MASS = 22.989770

    # Step 1: Define the sequences of the two peptides involved in the first disulfide bridge.
    # We hypothesize a typo in the C-terminal peptide ('GGEN' -> 'GEEN').
    peptide1_seq = "YDDMAACMK"
    peptide2_seq_hypothesized = "TQGCDEAEAGEGEEN"

    # Step 2: Calculate the neutral mass of each peptide.
    # Neutral Mass = (Sum of residue masses) + H2O
    mass_pep1 = sum(aa_masses[aa] for aa in peptide1_seq) + H2O_MASS
    mass_pep2 = sum(aa_masses[aa] for aa in peptide2_seq_hypothesized) + H2O_MASS

    # Step 3: Calculate the neutral mass of the disulfide-linked complex.
    # Mass_linked = Mass_pep1 + Mass_pep2 - 2*H
    linked_mass_neutral = mass_pep1 + mass_pep2 - (2 * H_MASS)

    # Step 4: Calculate the m/z for the hypothesized [M+H+Na]2+ adduct ion.
    # This ion has a total charge of +2.
    z = 2
    mz_value = (linked_mass_neutral + PROTON_MASS + NA_MASS) / z

    # Step 5: Print the results step-by-step
    print("--- Analysis of Disulfide Bridge 1 with Hypothesized Correction ---")
    print(f"Peptide 1 Sequence: {peptide1_seq}")
    print(f"Calculated neutral mass of Peptide 1: {mass_pep1:.3f}")
    print("\n")
    print(f"Hypothesized Peptide 2 Sequence (Corrected): {peptide2_seq_hypothesized}")
    print(f"Calculated neutral mass of Peptide 2: {mass_pep2:.3f}")
    print("\n")
    print(f"Neutral mass of linked peptides ({peptide1_seq} + {peptide2_seq_hypothesized} - 2H):")
    # To be fully transparent, we show the numbers in the equation
    print(f"{mass_pep1:.3f} + {mass_pep2:.3f} - {2 * H_MASS:.3f} = {linked_mass_neutral:.3f}")
    print("\n")
    print(f"Calculating m/z for the [M+H+Na]2+ adduct ion (charge z=2):")
    # Show numbers in the final equation
    print(f"({linked_mass_neutral:.3f} + {PROTON_MASS:.3f} + {NA_MASS:.3f}) / {z} = {mz_value:.3f}")
    print("\n")
    print(f"Final calculated m/z value: {mz_value:.3f}")


calculate_mz()