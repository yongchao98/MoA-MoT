import collections

def calculate_mz():
    """
    Calculates the m/z of a specific disulfide-linked peptide complex
    expected from the tryptic digest of active protein XER22.
    """

    # Monoisotopic masses of amino acid residues (mass loss of H2O is accounted for later)
    AA_MASSES = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }

    # Masses of other components
    H2O_MASS = 18.01056
    PROTON_MASS = 1.007276
    H_ATOM_MASS = 1.007825
    
    # Common modifications
    DEAMIDATION_MASS_SHIFT = 0.9840 # Mass difference for Q -> E or N -> D

    # Define the two peptides forming the disulfide bridge complex to be calculated.
    # This scenario assumes one missed cleavage for Peptide 1 and a deamidation of its glutamine (Q).
    peptide1_seq = "AKLAEQAERYDDMAACMK" # Contains Cys-23, with one missed cleavage
    peptide2_seq = "TQGCDEAEAGEGGEN"    # Contains Cys-153

    # --- Step 1: Calculate the mass of the first peptide ---
    peptide1_counts = collections.Counter(peptide1_seq)
    peptide1_mass = sum(AA_MASSES[aa] * count for aa, count in peptide1_counts.items())
    peptide1_mass += H2O_MASS # Add mass of terminal H and OH
    
    # Add mass for deamidation of the Glutamine (Q) residue
    peptide1_mass += DEAMIDATION_MASS_SHIFT

    # --- Step 2: Calculate the mass of the second peptide ---
    peptide2_counts = collections.Counter(peptide2_seq)
    peptide2_mass = sum(AA_MASSES[aa] * count for aa, count in peptide2_counts.items())
    peptide2_mass += H2O_MASS

    # --- Step 3: Calculate the mass of the disulfide-linked complex ---
    # The formation of a disulfide bridge removes two hydrogen atoms.
    disulfide_mass_loss = 2 * H_ATOM_MASS
    complex_mass = peptide1_mass + peptide2_mass - disulfide_mass_loss

    # --- Step 4: Calculate the m/z for a +3 charge state ---
    charge_state = 3
    mz_value = (complex_mass + charge_state * PROTON_MASS) / charge_state
    
    # --- Step 5: Print the results ---
    print(f"Analysis for Disulfide Bridge 1 with one missed cleavage and deamidation:")
    print("-" * 70)
    print(f"Peptide 1 Sequence: {peptide1_seq}")
    print(f"Calculated Mass of Peptide 1 (with deamidation): {round(peptide1_mass, 3)}")
    print(f"\nPeptide 2 Sequence: {peptide2_seq}")
    print(f"Calculated Mass of Peptide 2: {round(peptide2_mass, 3)}")
    print("-" * 70)
    print(f"Mass of Disulfide-Linked Complex [M]: {round(complex_mass, 3)}")
    print(f"Final m/z for [M + {charge_state}H]^{charge_state}+: {round(mz_value, 3)}")
    print("-" * 70)

calculate_mz()