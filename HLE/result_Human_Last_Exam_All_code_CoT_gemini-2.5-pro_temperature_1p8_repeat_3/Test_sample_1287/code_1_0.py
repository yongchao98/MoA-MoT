def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized, permethylated, and sodiated A2G2S2 glycan.

    The process involves:
    1. Calculating the mass of the initial A2G2S2 glycan.
    2. Calculating the mass change from converting two sialic acid carboxyl groups to amides.
    3. Calculating the mass added by permethylating all -OH and -NH sites.
    4. Summing these masses and adding the mass of a sodium ion.
    """

    # Monoisotopic atomic masses (Da)
    H = 1.007825
    C = 12.000000
    N = 14.003074
    O = 15.994915
    Na = 22.989770

    # Step 1: Calculate the mass of the unmodified, released A2G2S2 glycan
    # Composition: Hexose(5)HexNAc(4)Neu5Ac(2) -> 11 residues, 10 bonds
    # Formula = 5*C6H12O6 + 4*C8H15NO6 + 2*C11H19NO9 - 10*H2O
    # Final Formula: C84 H138 N6 O62
    num_C_base = 84
    num_H_base = 138
    num_N_base = 6
    num_O_base = 62
    mass_base_glycan = (num_C_base * C +
                        num_H_base * H +
                        num_N_base * N +
                        num_O_base * O)

    # Step 2: Calculate mass change from DMT-MM amidation
    # Reaction: R-COOH -> R-CONH2. Net change per group: -O + NH.
    # This occurs on two sialic acids.
    mass_change_amidation = 2 * (-O + N + H)

    # Step 3: Calculate mass change from permethylation
    # A methyl group (CH3) replaces a hydrogen (H). Net change: +CH2
    mass_per_methylation = C + 2 * H
    # Counting methylation sites on the amidated glycan:
    # A2G2 core structure: 23 sites (19 -OH, 4 -NH from 4 HexNAc)
    # Each modified sialic acid residue: 7 sites (4 -OH, 1 -NH from acetyl, 2 -NH from amide)
    num_sites_A2G2_part = 23
    num_sites_per_sialic_amide = 7
    total_methylation_sites = num_sites_A2G2_part + (2 * num_sites_per_sialic_amide)
    mass_change_methylation = total_methylation_sites * mass_per_methylation

    # Step 4: Calculate the final m/z for the [M+Na]+ ion
    mass_final_neutral = mass_base_glycan + mass_change_amidation + mass_change_methylation
    final_mz = mass_final_neutral + Na

    # Explanation
    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("The chemical modifications (amidation and permethylation) occur on functional groups")
    print("that are common to all three structures. The difference in sialic acid linkage")
    print("(alpha-2,3 vs. alpha-2,6) does not change the total number of atoms or methylation sites.")
    print("Therefore, all three isomers will have the identical mass-to-charge ratio (m/z).\n")

    # Print the detailed calculation
    print("Calculation of the m/z for the [M+Na]+ ion:")
    print("=" * 50)
    print(f"1. Mass of base glycan (C{num_C_base}H{num_H_base}N{num_N_base}O{num_O_base}): {mass_base_glycan:.4f} Da")
    print(f"2. Mass change from amidation of 2 sialic acids: {mass_change_amidation:.4f} Da")
    print(f"3. Mass change from permethylation of {total_methylation_sites} sites: {mass_change_methylation:.4f} Da")
    print(f"4. Mass of sodium adduct ([M+Na]+): {Na:.4f} Da")
    print("-" * 50)

    # Print the final equation
    print("Final Equation:")
    print(f"m/z = (Mass of Base Glycan) + (Mass Change Amidation) + (Mass Change Permethylation) + (Mass of Na)")
    print(f"m/z = {mass_base_glycan:.4f} + {mass_change_amidation:.4f} + {mass_change_methylation:.4f} + {Na:.4f}")
    print(f"m/z = {final_mz:.4f}\n")

    print("The expected mass for all three glycans is the same.")
    
    # Returning final value for the user request
    return f"{final_mz:.4f}"

final_answer = calculate_glycan_mass()
print(f"<<<{final_answer}>>>")