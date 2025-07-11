import math

def solve_disulfide_bridge_mz():
    """
    Calculates the m/z of a disulfide-linked peptide complex from protein XER22.
    """
    # Monoisotopic masses of amino acid residues (C5H9NO format)
    residue_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    # Mass constants
    H2O_MASS = 18.01056
    H_MASS = 1.007825
    PROTON_MASS = 1.007276

    # Peptides for the first disulfide bridge
    pep1_seq = 'YDDMAACMK'
    pep2_seq = 'TQGCDEAEAGEGGEN'

    print("Step 1: Calculate the mass of the first peptide.")
    mass_pep1 = sum(residue_masses[aa] for aa in pep1_seq) + H2O_MASS
    print(f"The sequence of peptide 1 is: {pep1_seq}")
    print(f"The calculated neutral mass of peptide 1 is: {mass_pep1:.5f} Da\n")

    print("Step 2: Calculate the mass of the second peptide and account for modification.")
    mass_pep2_original = sum(residue_masses[aa] for aa in pep2_seq) + H2O_MASS
    print(f"The sequence of peptide 2 is: {pep2_seq}")
    print(f"The original calculated neutral mass of peptide 2 is: {mass_pep2_original:.5f} Da")
    
    # A common modification is the deamidation of the C-terminal asparagine (N) to aspartic acid (D).
    # This often occurs during sample preparation or protein production.
    mass_n_residue = residue_masses['N']
    mass_d_residue = residue_masses['D']
    mass_pep2_deamidated = mass_pep2_original - mass_n_residue + mass_d_residue
    print("Considering a common deamidation of the C-terminal Asparagine (N) to Aspartic Acid (D):")
    print(f"The mass of deamidated peptide 2 is: {mass_pep2_original:.5f} - {mass_n_residue:.5f} + {mass_d_residue:.5f} = {mass_pep2_deamidated:.5f} Da\n")

    print("Step 3: Calculate the mass of the disulfide-linked peptide complex.")
    # The formation of a disulfide bond removes two hydrogen atoms.
    mass_linked_neutral = mass_pep1 + mass_pep2_deamidated - (2 * H_MASS)
    print("Mass = Mass(Pep1) + Mass(Pep2_deamidated) - 2 * Mass(H)")
    print(f"Mass = {mass_pep1:.5f} + {mass_pep2_deamidated:.5f} - 2 * {H_MASS:.5f} = {mass_linked_neutral:.5f} Da\n")

    print("Step 4: Calculate the m/z for the doubly charged ion (z=2).")
    charge = 2
    mz = (mass_linked_neutral + (charge * PROTON_MASS)) / charge
    print("m/z = (Neutral_Mass + charge * Mass(Proton)) / charge")
    # Outputting each number in the final equation as requested
    print(f"m/z = ({mass_linked_neutral:.5f} + {charge} * {PROTON_MASS:.5f}) / {charge}")
    print(f"The final calculated m/z value is: {round(mz, 3)}")

solve_disulfide_bridge_mz()