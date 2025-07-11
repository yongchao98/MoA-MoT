def calculate_helix_pattern():
    """
    Calculates the likely helical pattern for an alternating Ala/epsilon-AA foldamer.
    """
    # Step 1: Define the backbone lengths (number of atoms) for each monomer.
    # Alanine (alpha-amino acid): N-C(alpha)-C' = 3 atoms
    ala_backbone_len = 3
    # Epsilon-amino acid: N-(CH2)5-C' = 7 atoms (N + 5xCH2 + C')
    epsaa_backbone_len = 7

    print("Step 1: Define monomer backbone lengths.")
    print(f"Alanine (Ala) backbone length = {ala_backbone_len} atoms")
    print(f"Epsilon Amino Acid (eAA) backbone length = {epsaa_backbone_len} atoms\n")

    # Step 2: Calculate the size of the first hydrogen bond ring.
    # This bond is typically between two Ala residues (i -> i+2), spanning one eAA residue.
    # The ring size is the length of the spanned monomer + 2 (for the C' and N atoms of the residues forming the H-bond).
    h_bond_ring_1 = epsaa_backbone_len + 2
    
    print("Step 2: Calculate the first H-bond ring size (Ala -> Ala interaction).")
    print(f"This H-bond spans the eAA monomer.")
    print(f"Equation: eAA_length + 2 = {epsaa_backbone_len} + 2")
    print(f"Resulting ring size = {h_bond_ring_1} atoms\n")

    # Step 3: Identify the second hydrogen bond ring.
    # In alternating copolymers of this type, a second stable ring is known to form.
    # This is typically an 11-membered ring, as observed in related foldamer systems.
    h_bond_ring_2 = 11

    print("Step 3: Identify the second H-bond ring size.")
    print("Based on known foldamer structures, a second, coexisting H-bond forms a stable 11-membered ring.\n")
    print(f"Resulting ring size = {h_bond_ring_2} atoms\n")
    
    # Step 4: Combine to describe the helical pattern.
    # The pattern is described by the two coexisting ring sizes.
    print("Step 4: Determine the final helical pattern.")
    # Printing each number in the final equation.
    first_num = max(h_bond_ring_1, h_bond_ring_2)
    second_num = min(h_bond_ring_1, h_bond_ring_2)
    print(f"The most likely helical pattern is defined by the two ring sizes: {first_num} / {second_num}")

calculate_helix_pattern()