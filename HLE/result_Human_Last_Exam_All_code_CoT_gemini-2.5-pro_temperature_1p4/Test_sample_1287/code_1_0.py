import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


def calculate_mass():
    """
    Calculates the expected m/z for derivatized sialylated glycans.
    The script calculates the mass of the initial A2G2S2 glycan, then accounts for
    the mass changes from amidation and permethylation, and finally adds the
    mass of a sodium ion to find the m/z of the [M+Na]+ species.
    """
    # Monoisotopic atomic masses
    C_MASS = 12.000000
    H_MASS = 1.007825
    N_MASS = 14.003074
    O_MASS = 15.994915
    NA_MASS = 22.989770

    # Calculate masses of full monosaccharide molecules
    GLCNAC_MASS = 8 * C_MASS + 15 * H_MASS + 1 * N_MASS + 6 * O_MASS  # C8H15NO6
    MAN_MASS = 6 * C_MASS + 12 * H_MASS + 6 * O_MASS  # C6H12O6
    GAL_MASS = 6 * C_MASS + 12 * H_MASS + 6 * O_MASS  # C6H12O6
    NEU5AC_MASS = 11 * C_MASS + 19 * H_MASS + 1 * N_MASS + 9 * O_MASS # C11H19NO9
    WATER_MASS = 2 * H_MASS + 1 * O_MASS

    # --- Step 1: Calculate the mass of the starting glycan (A2G2S2) ---
    glycan_composition = {
        'GlcNAc': 4,
        'Man': 3,
        'Gal': 2,
        'Neu5Ac': 2
    }
    total_residues = sum(glycan_composition.values())
    num_bonds = total_residues - 1

    initial_glycan_mass = (glycan_composition['GlcNAc'] * GLCNAC_MASS +
                           glycan_composition['Man'] * MAN_MASS +
                           glycan_composition['Gal'] * GAL_MASS +
                           glycan_composition['Neu5Ac'] * NEU5AC_MASS -
                           num_bonds * WATER_MASS)

    # --- Step 2: Calculate the mass change from amidation ---
    # Reaction: -COOH -> -CONH2 per sialic acid. Net change: -O + N + H
    amidation_delta = -O_MASS + N_MASS + H_MASS
    num_amidations = glycan_composition['Neu5Ac']
    total_amidation_delta = num_amidations * amidation_delta

    # --- Step 3: Calculate the mass change from permethylation ---
    # Mass change per methylation (-H -> -CH3) is +CH2
    methylation_unit_mass = C_MASS + 2 * H_MASS

    # Count total methylation sites
    oh_per_residue = {'GlcNAc': 3, 'Man': 4, 'Gal': 4, 'Neu5Ac': 4}
    total_oh_groups = sum(glycan_composition[res] * oh_per_residue[res] for res in oh_per_residue)
    protons_per_amide = 2
    total_amide_protons = num_amidations * protons_per_amide
    total_methylation_sites = total_oh_groups + total_amide_protons
    total_methylation_delta = total_methylation_sites * methylation_unit_mass

    # --- Step 4: Calculate final neutral mass and m/z of the sodiated ion ---
    final_neutral_mass = initial_glycan_mass + total_amidation_delta + total_methylation_delta
    final_mz = final_neutral_mass + NA_MASS

    # --- Print the explanation and final result ---
    print("The three glycans you listed, A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2, are isomers.")
    print("This means they have the same chemical composition and therefore the same mass. The chemical modifications")
    print("(amidation and permethylation) also affect the same number of functional groups in each isomer.")
    print("As a result, they will all have the same final m/z value.")
    print("\nHere is the step-by-step calculation for the final m/z of the singly-sodiated ion [M+Na]+:\n")
    
    print(f"1. Mass of starting glycan (A2G2S2): {initial_glycan_mass:.4f} Da")
    
    print(f"\n2. Mass change from amidation of {num_amidations} sialic acids:")
    print(f"   - Change per reaction (-O+NH): {amidation_delta:.4f} Da")
    print(f"   - Total change: {num_amidations} * {amidation_delta:.4f} = {total_amidation_delta:.4f} Da")
    
    print(f"\n3. Mass change from permethylation:")
    print(f"   - Total methylation sites (-OH + -NH2 protons): {total_oh_groups} + {total_amide_protons} = {total_methylation_sites}")
    print(f"   - Change per methylation (+CH2): {methylation_unit_mass:.4f} Da")
    print(f"   - Total change: {total_methylation_sites} * {methylation_unit_mass:.4f} = {total_methylation_delta:.4f} Da")
    
    print(f"\n4. Mass of Sodium adduct (Na+): {NA_MASS:.4f} Da")

    print("\n-------------------------------------------------------------")
    print("Final m/z Calculation:")
    print("Mass(A2G2S2) + Change(Amidation) + Change(Permethylation) + Mass(Na+) = Final m/z")
    print(f"{initial_glycan_mass:.4f} Da + ({total_amidation_delta:.4f} Da) + {total_methylation_delta:.4f} Da + {NA_MASS:.4f} Da = {final_mz:.4f} m/z")
    print("-------------------------------------------------------------")

    return final_mz

# Execute the function and capture the output
final_mass_to_charge = calculate_mass()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

# Print the final answer in the specified format
print(f"<<<{final_mass_to_charge:.4f}>>>")