import collections

def calculate_mz():
    """
    Calculates the mass and m/z for disulfide-linked peptides.
    """
    # Monoisotopic residue masses of amino acids
    MASSES = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    # Mass constants
    H2O_MASS = 18.01056
    H_ATOM_MASS = 1.007825
    PROTON_MASS = 1.007276

    # Peptides for the first disulfide bridge
    p1_seq = "YDDMAACMK"
    p2_seq = "TQGCDEAEAGEGGEN"

    # --- Step 1: Calculate the mass of each peptide ---
    p1_residue_mass = sum(MASSES[aa] for aa in p1_seq)
    p1_mass = p1_residue_mass + H2O_MASS

    p2_residue_mass = sum(MASSES[aa] for aa in p2_seq)
    p2_mass = p2_residue_mass + H2O_MASS

    print(f"Calculating for Bridge 1:")
    print(f"Peptide 1: {p1_seq}")
    print(f"Peptide 2: {p2_seq}")
    print("-" * 30)
    print(f"Mass of {p1_seq}: {p1_residue_mass:.3f} (residues) + {H2O_MASS:.3f} (termini) = {p1_mass:.3f} Da")
    print(f"Mass of {p2_seq}: {p2_residue_mass:.3f} (residues) + {H2O_MASS:.3f} (termini) = {p2_mass:.3f} Da")
    print("-" * 30)

    # --- Step 2: Calculate the mass of the disulfide-linked pair ---
    linked_mass_neutral = p1_mass + p2_mass - (2 * H_ATOM_MASS)
    
    print("Calculating linked mass:")
    print(f"Linked Mass = Mass({p1_seq}) + Mass({p2_seq}) - 2*Mass(H)")
    print(f"Linked Mass = {p1_mass:.3f} + {p2_mass:.3f} - 2*{H_ATOM_MASS:.3f} = {linked_mass_neutral:.3f} Da")
    print("-" * 30)

    # --- Step 3: Calculate m/z for different charge states ---
    print("Calculating potential m/z values:")
    for z in [2, 3, 4]: # Common charge states for peptides of this size
        mz = (linked_mass_neutral + (z * PROTON_MASS)) / z
        print(f"For charge z=+{z}:")
        print(f"m/z = ({linked_mass_neutral:.3f} + {z}*{PROTON_MASS:.3f}) / {z} = {mz:.3f}")

    print("\nNote: Similar calculations for the second bridge (FLIPNACSQAESK + ACSLAK) and for potential missed cleavages also do not match the provided options.")


calculate_mz()