def calculate_glycan_mass():
    """
    Calculates the theoretical monoisotopic mass of a specific permethylated and
    amidated glycan observed as a sodiated ion in LC-MS.
    """

    # Monoisotopic atomic masses (Da)
    C_mass = 12.000000
    H_mass = 1.007825
    O_mass = 15.994915
    N_mass = 14.003074
    Na_mass = 22.989770

    # Step 1: Calculate the mass of the native A2G2S2 glycan.
    # Composition: 3 Man, 4 GlcNAc, 2 Gal, 2 Neu5Ac
    # Formula of monosaccharides:
    # Man/Gal (Hexose): C6H12O6
    # GlcNAc (HexNAc): C8H15NO6
    # Neu5Ac (Sialic Acid): C11H19NO9
    # Total residues (n) = 3 + 4 + 2 + 2 = 11
    # To form the polysaccharide, (n-1) = 10 molecules of H2O are removed.

    num_C = (3 * 6) + (4 * 8) + (2 * 6) + (2 * 11)
    num_H = (3 * 12) + (4 * 15) + (2 * 12) + (2 * 19) - (10 * 2)
    num_N = (4 * 1) + (2 * 1)
    num_O = (3 * 6) + (4 * 6) + (2 * 6) + (2 * 9) - (10 * 1)

    print(f"The chemical formula for the native glycan is C{num_C}H{num_H}N{num_N}O{num_O}.")

    native_mass = (num_C * C_mass) + (num_H * H_mass) + (num_N * N_mass) + (num_O * O_mass)
    print("Calculating mass of the native glycan...")
    print(f"Mass = ({num_C} * {C_mass:.6f}) + ({num_H} * {H_mass:.6f}) + ({num_N} * {N_mass:.6f}) + ({num_O} * {O_mass:.6f}) = {native_mass:.4f} Da")
    print("-" * 30)

    # Step 2: Apply amidation of the two sialic acids.
    # The reaction -COOH -> -CONH2 replaces one Oxygen with an NH group.
    num_amidations = 2
    amidation_mass_change = N_mass + H_mass - O_mass
    amidated_mass = native_mass + (num_amidations * amidation_mass_change)
    print("Step 1: Amidation (2x)")
    print(f"Mass change per amidation = mass(N) + mass(H) - mass(O) = {amidation_mass_change:.4f} Da")
    print(f"Mass after amidation = {native_mass:.4f} + (2 * {amidation_mass_change:.4f}) = {amidated_mass:.4f} Da")
    print("-" * 30)

    # Step 3: Apply permethylation.
    # This adds a CH2 group (replaces H with CH3) for each available site.
    # For this glycan structure, there are 32 methylation sites.
    num_methylations = 32
    methylation_mass_change = C_mass + (2 * H_mass)
    permethylated_mass = amidated_mass + (num_methylations * methylation_mass_change)
    print(f"Step 2: Permethylation ({num_methylations}x)")
    print(f"Mass change per methylation = mass(C) + 2*mass(H) = {methylation_mass_change:.4f} Da")
    print(f"Mass after permethylation = {amidated_mass:.4f} + ({num_methylations} * {methylation_mass_change:.4f}) = {permethylated_mass:.4f} Da")
    print("-" * 30)

    # Step 4: Add mass of sodium for the singly-sodiated ion [M+Na]+.
    final_mass = permethylated_mass + Na_mass
    print("Step 3: Sodiated Ion Formation [M+Na]+")
    print(f"Final observed mass = {permethylated_mass:.4f} + mass(Na) = {final_mass:.4f} Da")
    print("-" * 30)
    
    print("\nSince all three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers, they will all have the same final mass.")
    print(f"The expected mass for the [M+Na]+ ion for all three glycans is approximately {final_mass:.3f} Da.")

if __name__ == "__main__":
    calculate_glycan_mass()