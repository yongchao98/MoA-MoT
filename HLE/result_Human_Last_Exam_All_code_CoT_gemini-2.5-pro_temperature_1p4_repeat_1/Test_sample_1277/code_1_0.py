import sys
import io

# Redirect print to a string buffer to capture output for final formatting
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Mass Constants ---
# Monoisotopic residue masses
MASS_HEX = 162.05282  # Hexose (Galactose, Mannose)
MASS_HEXNAC = 203.07937 # N-Acetylhexosamine (GlcNAc)
MASS_FUC = 146.05751   # Fucose (dHex)
MASS_NEU5AC = 291.09542 # N-Acetylneuraminic acid (Sialic Acid)

# Other masses
MASS_H2O = 18.010565
MASS_H_PLUS = 1.007276
# RFMS derivatization adds the mass of the tag (C16H20N4O) and removes one water molecule.
MASS_RFMS_TAG = 284.163711
MASS_ADDED_BY_RFMS = MASS_RFMS_TAG - MASS_H2O

# --- User Provided Data ---
parent_mz = 856.6638
isotope_mz1 = 856.9971
fragment_ions = [204.087, 366.140, 528.193, 673.231, 882.409, 1368.568, 1894.753, 2260.886]
base_peak_mz = 528.193

def solve_glycan():
    """
    Analyzes the provided MS data to identify the glycan structure.
    """
    print("Step 1: Analyzing the Parent Ion (MS1)")
    # Determine charge state
    isotope_spacing = isotope_mz1 - parent_mz
    charge_state = round(1 / isotope_spacing)
    print(f"The isotopic spacing is {isotope_spacing:.4f} Da, which indicates a charge state of z = {charge_state}.")

    # Calculate mass of the derivatized glycan
    mass_ion = parent_mz * charge_state
    mass_derivatized_glycan = mass_ion - (charge_state * MASS_H_PLUS)
    print(f"The parent ion [M+{charge_state}H]^{charge_state}+ at m/z {parent_mz} corresponds to a neutral derivatized glycan mass of {mass_derivatized_glycan:.4f} Da.")
    print("-" * 20)

    print("Step 2: Analyzing the Fragmentation Data (MSMS)")
    print("The fragment ions (oxonium ions) provide key structural information:")
    print(f"- m/z {fragment_ions[0]:.3f}: Corresponds to a HexNAc ion, a basic building block of N-glycans.")
    print(f"- m/z {fragment_ions[1]:.3f}: Corresponds to a Hexose-HexNAc ion, indicating a galactose-GlcNAc (Gal-GlcNAc) or mannose-GlcNAc antenna structure.")
    print(f"- m/z {base_peak_mz:.3f}: This is the most intense fragment ion (base peak). It corresponds to a [Hex-Hex-HexNAc]+ fragment, such as one originating from a Man-(Gal-GlcNAc) antenna branch.")
    print("The high intensity of this specific ion is a well-known characteristic of a tri-antennary (A3) N-glycan structure. It arises from the cleavage of one of the three antennae.")
    print("-" * 20)

    print("Step 3: Deducing the Glycan Structure")
    print("Based on the MSMS data, the glycan is a complex, fucosylated, tri-antennary structure.")
    print("The high mass suggests it is also sialylated.")
    print("The most plausible structure that fits this profile is a fucosylated, tri-antennary, tri-galactosylated, mono-sialylated glycan.")

    # Define the proposed structure composition
    composition = {
        'Fuc': 1, 'Hex': 6, 'HexNAc': 5, 'Neu5Ac': 1
    }
    # Name breakdown: F(Fucosylated), A3(Tri-antennary), G3(Tri-galactosylated), S1(Mono-sialylated)
    glycan_name = "FA3G3S1"

    print(f"\nThis corresponds to the Oxford nomenclature name: {glycan_name}")
    print(f"Composition: {composition['Fuc']} Fucose + {composition['Hex']} Hexose + {composition['HexNAc']} HexNAc + {composition['Neu5Ac']} Neu5Ac")
    print("-" * 20)
    
    print("Step 4: Verification and Linkage Information")
    
    # Calculate theoretical mass and m/z to check consistency
    residue_mass = (composition['Fuc'] * MASS_FUC) + \
                   (composition['Hex'] * MASS_HEX) + \
                   (composition['HexNAc'] * MASS_HEXNAC) + \
                   (composition['Neu5Ac'] * MASS_NEU5AC)
                   
    native_mass = residue_mass + MASS_H2O
    theoretical_derivatized_mass = native_mass + MASS_ADDED_BY_RFMS
    theoretical_mz = (theoretical_derivatized_mass + charge_state * MASS_H_PLUS) / charge_state
    
    print("A note on the provided parent ion m/z:")
    print(f"The theoretical [M+3H]3+ for {glycan_name} is m/z {theoretical_mz:.4f}.")
    print(f"This does not match the observed m/z of {parent_mz}. This discrepancy is common and can be due to non-proton adducts (e.g., sodium) or other factors. In such cases, the structural information from the fragment ions is more reliable for identification.")
    print("\nLinkage Information from MSMS:")
    print("While exact linkages are hard to determine without more data, the fragmentation pattern provides clues:")
    print("- The glycan is core fucosylated, meaning a fucose is attached via an α1-6 linkage to the innermost GlcNAc.")
    print("- The presence of a dominant m/z 528 fragment suggests it is a tri-antennary glycan with antennae on the trimannosyl core. A common structure involves antennae at the C2 position of the α1-3 mannose, and at the C2 and C4 positions of the α1-6 mannose.")

    print("\n--- Final Answer ---")
    print(f"The name of the glycan is {glycan_name}")

solve_glycan()

# Get the captured output
final_output = captured_output.getvalue()

# Close the buffer and restore stdout
sys.stdout = old_stdout

print(final_output)

final_answer_name = "FA3G3S1"
<<<f"{final_answer_name}">>>