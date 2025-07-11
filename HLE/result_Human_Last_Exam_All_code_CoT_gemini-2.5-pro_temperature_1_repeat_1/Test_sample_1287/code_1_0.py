import sys

def solve():
    """
    Calculates the expected m/z for modified sialylated glycans.
    """
    # Define monoisotopic masses of elements and sodium
    mass_H = 1.007825
    mass_C = 12.000000
    mass_N = 14.003074
    mass_O = 15.994915
    mass_Na = 22.989770

    # Define formulas for glycan residues (mass after loss of H2O)
    # Hexose (Galactose, Mannose) residue: C6H10O5
    hex_residue_mass = 6 * mass_C + 10 * mass_H + 5 * mass_O
    # N-acetylglucosamine (GlcNAc) residue: C8H13NO5
    glcnac_residue_mass = 8 * mass_C + 13 * mass_H + 1 * mass_N + 5 * mass_O
    # N-acetylneuraminic acid (Sialic acid, Neu5Ac) residue: C11H17NO8
    neu5ac_residue_mass = 11 * mass_C + 17 * mass_H + 1 * mass_N + 8 * mass_O

    # Mass of water to add for the non-reducing end hydroxyls
    mass_H2O = 2 * mass_H + mass_O

    # Step 1: Calculate the mass of the native A2G2S2 glycan
    # Composition: 3 Man + 2 Gal + 4 GlcNAc + 2 Neu5Ac
    # This simplifies to 5 Hexose + 4 GlcNAc + 2 Neu5Ac
    num_hex = 5
    num_glcnac = 4
    num_neu5ac = 2

    native_glycan_mass = (num_hex * hex_residue_mass +
                          num_glcnac * glcnac_residue_mass +
                          num_neu5ac * neu5ac_residue_mass +
                          mass_H2O)

    # Step 2: Calculate mass change from amidation
    # The reaction converts -COOH on each sialic acid to -CONH2.
    # This corresponds to replacing an -OH group with an -NH2 group.
    # Net atomic change is -O, +N, +H per sialic acid.
    num_sialic_acids = 2
    amidation_mass_change = num_sialic_acids * (mass_N + mass_H - mass_O)
    mass_after_amidation = native_glycan_mass + amidation_mass_change

    # Step 3: Calculate mass change from permethylation
    # Permethylation replaces H on all OH and NH groups with a methyl group (CH3).
    # The mass change per site = Mass(CH3) - Mass(H) = C + 2H
    methylation_unit_mass_gain = mass_C + 2 * mass_H

    # Count the total number of methylation sites on the amidated glycan.
    # The count is based on the standard A2G2S2 N-glycan structure.
    # Reducing end GlcNAc: 3 OH + 1 NH = 4 sites
    # Core GlcNAc: 2 OH + 1 NH = 3 sites
    # Branching Man: 2 OH = 2 sites
    # Antennae Man (x2): 3 OH each = 6 sites
    # Antennae GlcNAc (x2): 2 OH + 1 NH each = 6 sites
    # Antennae Gal (x2): 3 OH each = 6 sites
    # Amidated Sia (x2): 4 OH + 1 NH (acetyl) + 2 NH (amide) each = 14 sites
    # Total sites = 4 + 3 + 2 + 6 + 6 + 6 + 14 = 41 sites
    num_methylation_sites = 41
    permethylation_mass_gain = num_methylation_sites * methylation_unit_mass_gain

    # Step 4: Calculate the final mass of the modified, permethylated glycan
    final_neutral_mass = mass_after_amidation + permethylation_mass_gain

    # Step 5: Calculate m/z for the singly-sodiated ion [M+Na]+
    final_mz = final_neutral_mass + mass_Na

    # Step 6: Print the results
    # The final mass is identical for all three linkage isomers.
    glycans = ["A2G(4)2S(3)2", "A2G(4)S(3)S(6)", "A2G(4)2S(6)2"]
    
    print("The chemical modifications and resulting mass are independent of the sialic acid linkage (alpha-2,3 vs. alpha-2,6).")
    print("Therefore, all three isomers will have the same theoretical m/z.")
    print("-" * 60)
    
    # Use a variable to store the final result for the "answer" block
    final_answer_value = 0.0

    for i, glycan in enumerate(glycans):
        print(f"The expected singly-sodiated m/z for {glycan} is calculated as:")
        # The prompt requires printing each number in the final equation.
        print(f"  m/z = Mass(Native Glycan) + Δ(Amidation) + Δ(Permethylation) + Mass(Na+)")
        print(f"  m/z = {native_glycan_mass:.4f} + ({amidation_mass_change:.4f}) + {permethylation_mass_gain:.4f} + {mass_Na:.4f}")
        print(f"  m/z = {final_mz:.4f}")
        print()
        if i == len(glycans) - 1:
            final_answer_value = final_mz

    # The final answer is the same for all three.
    sys.stdout.write(f'<<<{final_answer_value:.4f}>>>')

solve()