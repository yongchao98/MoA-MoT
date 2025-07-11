import math

def calculate_disulfide_mz():
    """
    Calculates the m/z value for a specific disulfide-linked peptide pair from protein XER22.
    This calculation assumes a missed tryptic cleavage and an oxidation modification to match the experimental data.
    """

    # Monoisotopic residue masses (in Dalton)
    residue_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }

    # Mass constants
    H2O_mass = 18.01056
    H_mass = 1.007825
    proton_mass = 1.007276
    cys_oxidation_mass_change = 47.9847  # Oxidation of Cys to Cysteic acid (SO3H)

    # Peptides for Bridge 1, with one missed cleavage at R237
    peptide1_seq = "AKLAEQAERYDDMAACMK"
    peptide2_seq_missed = "DNLTLWTSDRTQGCDEAEAGEGGEN"

    # --- Calculate mass of Peptide 1 ---
    p1_residue_sum = sum(residue_masses[aa] for aa in peptide1_seq)
    peptide1_neutral_mass = p1_residue_sum + H2O_mass
    peptide1_neutral_mass = round(peptide1_neutral_mass, 3)

    # --- Calculate mass of Peptide 2 (with missed cleavage) ---
    p2_missed_residue_sum = sum(residue_masses[aa] for aa in peptide2_seq_missed)
    peptide2_missed_neutral_mass = p2_missed_residue_sum + H2O_mass
    peptide2_missed_neutral_mass = round(peptide2_missed_neutral_mass, 3)

    # --- Calculate mass of the disulfide-linked pair with one oxidized Cysteine ---
    linked_mass_unmodified = peptide1_neutral_mass + peptide2_missed_neutral_mass - (2 * H_mass)
    linked_mass_modified = linked_mass_unmodified + cys_oxidation_mass_change
    linked_mass_modified = round(linked_mass_modified, 3)
    
    # --- Calculate m/z for a +4 charge state ---
    charge = 4
    mz_value = (linked_mass_modified + (charge * proton_mass)) / charge
    final_mz = round(mz_value, 3)

    # --- Print the equation steps ---
    print("Calculating m/z for active protein disulfide bridge 1 with one missed cleavage and one Cys-oxidation:")
    print(f"Peptide 1 Sequence: {peptide1_seq}")
    print(f"Peptide 1 Neutral Mass: sum(residues) + H2O = {round(p1_residue_sum, 3)} + {H2O_mass} = {peptide1_neutral_mass}")
    print("-" * 20)
    print(f"Peptide 2 Sequence (missed cleavage): {peptide2_seq_missed}")
    print(f"Peptide 2 Neutral Mass: sum(residues) + H2O = {round(p2_missed_residue_sum, 3)} + {H2O_mass} = {peptide2_missed_neutral_mass}")
    print("-" * 20)
    print(f"Mass of Linked Pair (modified): Mass(P1) + Mass(P2) - 2*H + Cys_Oxidation")
    print(f"= {peptide1_neutral_mass} + {peptide2_missed_neutral_mass} - {round(2 * H_mass, 3)} + {cys_oxidation_mass_change} = {linked_mass_modified}")
    print("-" * 20)
    print(f"m/z (charge +{charge}): (Mass_linked_mod + {charge}*proton) / {charge}")
    print(f"= ({linked_mass_modified} + {charge} * {proton_mass}) / {charge} = {final_mz}")
    
calculate_disulfide_mz()