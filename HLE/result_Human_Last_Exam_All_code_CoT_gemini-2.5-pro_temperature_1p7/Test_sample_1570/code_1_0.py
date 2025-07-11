import sys

def solve():
    """
    This function calculates the m/z of a disulfide-linked peptide complex from Protein XER22,
    assuming one missed cleavage, one phosphorylation event, and a +3 charge state.
    """
    # Monoisotopic masses of amino acid residues (AA - H2O)
    aa_mass = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }

    # Other relevant masses
    mass_h2o = 18.01056        # Mass of water
    mass_h = 1.007825          # Mass of a hydrogen atom
    mass_proton = 1.007276     # Mass of a proton (for m/z calculation)
    mass_phosphate = 79.96633  # Mass of HPO3 (phosphorylation)

    # --- Step 1: Define peptide sequences for Bridge 1 ---
    # Peptide 1 from a tryptic digest with one missed cleavage, contains Cys from 'MAACM'
    p1_seq = "AKLAEQAERYDDMAACMK"
    # Peptide 2 from a tryptic digest, contains Cys from 'TQGCDEAEAGEG'
    p2_seq = "TQGCDEAEAGEGGEN"

    # --- Step 2: Calculate the mass of each peptide ---
    # Calculate mass by summing residue masses and adding the mass of one water molecule.
    p1_residue_mass = sum(aa_mass[aa] for aa in p1_seq)
    p1_mass = p1_residue_mass + mass_h2o
    p1_mass = round(p1_mass, 3)

    p2_residue_mass = sum(aa_mass[aa] for aa in p2_seq)
    p2_mass = p2_residue_mass + mass_h2o
    p2_mass = round(p2_mass, 3)
    
    # --- Step 3: Account for phosphorylation ---
    # The peptide from HEK293F cells may be phosphorylated. P1 contains a 'Y' residue.
    p1_mass_phosphorylated = p1_mass + mass_phosphate
    p1_mass_phosphorylated = round(p1_mass_phosphorylated, 3)
    
    # --- Step 4: Calculate mass of the disulfide-linked complex ---
    # The mass is the sum of the two peptides minus two hydrogens.
    linked_mass_neutral = p1_mass_phosphorylated + p2_mass - (2 * mass_h)
    linked_mass_neutral = round(linked_mass_neutral, 3)

    # --- Step 5: Calculate the m/z value for a +3 charge state ---
    charge = 3
    mz_value = (linked_mass_neutral + (charge * mass_proton)) / charge
    mz_value = round(mz_value, 3)
    
    # --- Output the results ---
    print(f"Analysis for Disulfide Bridge 1 with Phosphorylation")
    print(f"--------------------------------------------------")
    print(f"Peptide 1 Sequence: {p1_seq}")
    print(f"Peptide 2 Sequence: {p2_seq}")
    print(f"Mass of Peptide 1 (unmodified): {p1_mass} Da")
    print(f"Mass of Peptide 2: {p2_mass} Da")
    print(f"Mass of phosphorylated Peptide 1: {p1_mass} + {round(mass_phosphate, 3)} = {p1_mass_phosphorylated} Da")
    print(f"Neutral mass of linked complex: {p1_mass_phosphorylated} + {p2_mass} - {round(2 * mass_h, 3)} = {linked_mass_neutral} Da")
    print(f"Final m/z for a +3 ion: ({linked_mass_neutral} + {charge} * {round(mass_proton, 3)}) / {charge} = {mz_value}")

solve()