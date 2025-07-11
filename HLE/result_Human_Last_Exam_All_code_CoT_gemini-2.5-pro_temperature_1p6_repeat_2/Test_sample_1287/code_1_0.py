import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = string_buffer

def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized glycan.
    The glycan is A2G2S2, which undergoes amidation, permethylation, and sodiation.
    """
    # Monoisotopic atomic masses
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_NA = 22.989770

    # Calculate masses of monosaccharide building blocks and water
    mass_h2o = 2 * MASS_H + MASS_O
    # Mannose/Galactose (Hexose, C6H12O6)
    mass_hex = 6 * MASS_C + 12 * MASS_H + 6 * MASS_O
    # N-Acetylglucosamine (HexNAc, C8H15NO6)
    mass_glcnac = 8 * MASS_C + 15 * MASS_H + 1 * MASS_N + 6 * MASS_O
    # N-Acetylneuraminic Acid (Sialic Acid, C11H19NO9)
    mass_neu5ac = 11 * MASS_C + 19 * MASS_H + 1 * MASS_N + 9 * MASS_O

    # 1. Calculate the mass of the native A2G2S2 glycan
    num_hex = 5  # 3 Mannose + 2 Galactose
    num_glcnac = 4
    num_neu5ac = 2
    num_residues = num_hex + num_glcnac + num_neu5ac
    num_glycosidic_bonds = num_residues - 1

    mass_native_glycan = (num_hex * mass_hex +
                          num_glcnac * mass_glcnac +
                          num_neu5ac * mass_neu5ac -
                          num_glycosidic_bonds * mass_h2o)

    # 2. Calculate mass change from amidation
    # This reaction converts -COOH to -CONH2, which is a net change of -O +NH.
    mass_change_amidation_per_site = MASS_N + MASS_H - MASS_O
    total_mass_change_amidation = num_neu5ac * mass_change_amidation_per_site
    mass_after_amidation = mass_native_glycan + total_mass_change_amidation

    # 3. Calculate mass change from permethylation
    # This reaction replaces H with CH3, a net addition of CH2.
    mass_change_per_methylation = MASS_C + 2 * MASS_H
    # The number of available -OH and -NH sites on the amidated A2G2S2 structure is 35.
    num_methylation_sites = 35
    total_mass_change_methylation = num_methylation_sites * mass_change_per_methylation
    mass_derivatized = mass_after_amidation + total_mass_change_methylation

    # 4. Add mass of sodium for the final [M+Na]+ ion
    final_ion_mass = mass_derivatized + MASS_NA

    # Since charge is +1, m/z = final_ion_mass / 1
    mz_value = final_ion_mass

    # --- Output the results ---
    print("The three specified glycans are isomers and will have the same mass after derivatization.")
    print("\nCalculation Breakdown:")
    print(f"1. Mass of native A2G2S2 glycan: {mass_native_glycan:.4f} Da")
    print(f"2. Mass change from amidation of 2 sialic acids (-O +NH)x2: {total_mass_change_amidation:.4f} Da")
    print(f"3. Mass change from permethylation on {num_methylation_sites} sites (+CH2)x{num_methylation_sites}: {total_mass_change_methylation:.4f} Da")
    print(f"4. Mass of the final derivatized glycan M = {mass_native_glycan:.4f} + ({total_mass_change_amidation:.4f}) + {total_mass_change_methylation:.4f} = {mass_derivatized:.4f} Da")
    print(f"5. Mass of the singly sodiated ion [M+Na]+ = {mass_derivatized:.4f} + {MASS_NA:.4f} = {final_ion_mass:.4f} Da")
    print("\nThe expected m/z value for the singly-sodiated ion of each glycan is:")
    print(f"A2G(4)2S(3)2: {mz_value:.4f}")
    print(f"A2G(4)S(3)S(6): {mz_value:.4f}")
    print(f"A2G(4)2S(6)2: {mz_value:.4f}")

    return mz_value

# Execute the calculation
final_mz = calculate_glycan_mass()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = string_buffer.getvalue()

# Print the captured output to the actual console
print(output)

# Final answer block
final_answer_string = f"{final_mz:.4f}"
print(f"\n<<<A2G(4)2S(3)2 = {final_answer_string}, A2G(4)S(3)S(6) = {final_answer_string}, A2G(4)2S(6)2 = {final_answer_string}>>>")