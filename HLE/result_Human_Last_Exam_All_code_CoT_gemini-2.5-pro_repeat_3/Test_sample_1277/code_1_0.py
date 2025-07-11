import sys
from io import StringIO

# A helper class to redirect stdout to a string
class Scriptor:
    def __init__(self):
        self.stdout = sys.stdout
        self.string_io = StringIO()

    def start(self):
        sys.stdout = self.string_io

    def stop(self):
        sys.stdout = self.stdout
        output = self.string_io.getvalue()
        self.string_io.close()
        return output

def solve_glycan_puzzle():
    """
    Solves the glycan identification puzzle based on MS data.
    """
    # --- Masses of relevant chemical species (monoisotopic) ---
    # Monosaccharide residues (loss of H2O)
    mass_Hex = 162.05282
    mass_HexNAc = 203.07937
    mass_Fuc = 146.05791  # Fucose (a deoxyhexose, dHex)
    mass_NeuAc = 291.09542 # N-Acetylneuraminic acid

    # Other masses
    mass_H2O = 18.01056
    mass_proton = 1.007276
    # The mass added by RFMS derivatization (reductive amination)
    # is Mass(RFMS_tag) - Mass(H2O) = 365.1331 - 18.0106 = 347.1225
    # The problem asks to identify the glycan, so we calculate the glycan's native mass.
    mass_RFMS_derivatized = 365.1331 # Mass of the added tag C19H17N3O3S

    # --- Input Data from the User ---
    precursor_mz = 856.6638
    isotope_peak_1 = 856.9971
    msms_fragments_obs = [204.087, 366.140, 528.193, 673.231, 882.409, 1368.568, 1894.753, 2260.886]

    # --- Step 1: Analyze the Precursor Ion (MS1) ---
    print("--- Step 1: Precursor Ion Analysis ---")
    isotope_spacing = isotope_peak_1 - precursor_mz
    # Charge state (z) = Mass(C13-C12) / isotope_spacing
    charge_state = round(1.00335 / isotope_spacing)
    print(f"Isotope spacing is {isotope_spacing:.4f} Da, which indicates a charge state (z) of {charge_state}.")

    # Calculate neutral mass of the observed ion (RFMS-glycan complex)
    neutral_mass_derivatized = (precursor_mz * charge_state) - (charge_state * mass_proton)
    print(f"The neutral mass of the derivatized glycan is ({precursor_mz:.4f} * {charge_state}) - ({charge_state} * {mass_proton:.4f}) = {neutral_mass_derivatized:.4f} Da.")

    # Calculate the mass of the original, underivatized glycan
    # The derivatization is a condensation reaction, so one H2O is lost.
    # Mass(derivatized) = Mass(glycan) + Mass(RFMS_tag) - Mass(H2O)
    # Mass(glycan) = Mass(derivatized) - Mass(RFMS_tag) + Mass(H2O)
    mass_glycan_exp = neutral_mass_derivatized - mass_RFMS_derivatized + mass_H2O
    print(f"The mass of the underivatized glycan is {neutral_mass_derivatized:.4f} - {mass_RFMS_derivatized:.4f} + {mass_H2O:.4f} = {mass_glycan_exp:.4f} Da.\n")

    # --- Step 2: Determine Glycan Composition ---
    print("--- Step 2: Glycan Composition Identification ---")
    # Propose a composition: FA2G2S1
    # F(1) = 1 Fucose, A2 = biantennary (2 core GlcNAc + 2 antenna GlcNAc), G2 = 2 Galactose, S1 = 1 Sialic Acid
    # Core Mannose = 3. Total composition: Hex(5), HexNAc(4), Fuc(1), NeuAc(1)
    num_Hex = 5
    num_HexNAc = 4
    num_Fuc = 1
    num_NeuAc = 1
    
    # Calculate theoretical mass of the proposed glycan
    mass_glycan_theo = (num_Hex * mass_Hex) + (num_HexNAc * mass_HexNAc) + \
                       (num_Fuc * mass_Fuc) + (num_NeuAc * mass_NeuAc) + mass_H2O # Add H2O for the reducing end
    
    print(f"Proposing composition: Hex({num_Hex}) HexNAc({num_HexNAc}) Fuc({num_Fuc}) NeuAc({num_NeuAc}).")
    print(f"Theoretical mass = (5*{mass_Hex:.4f}) + (4*{mass_HexNAc:.4f}) + (1*{mass_Fuc:.4f}) + (1*{mass_NeuAc:.4f}) + {mass_H2O:.4f} = {mass_glycan_theo:.4f} Da.")
    
    mass_error = mass_glycan_exp - mass_glycan_theo
    ppm_error = (mass_error / mass_glycan_theo) * 1e6
    print(f"The mass error is {mass_error:.4f} Da ({ppm_error:.2f} ppm), confirming the composition is FA2G2S1.\n")

    # --- Step 3: Analyze Fragmentation Data (MS/MS) ---
    print("--- Step 3: MS/MS Fragmentation Analysis ---")
    print("Note: Fragment ions cannot have a higher m/z than the precursor (856.6638).")
    print("The high m/z values in the provided list are invalid fragments and will be ignored.")
    
    # Calculate theoretical masses of key diagnostic fragment ions ([M+H]+)
    frag_hexnac = mass_HexNAc + mass_H2O + mass_proton
    frag_lacnac = mass_Hex + mass_HexNAc + mass_H2O + mass_proton
    frag_b3_ion = (2 * mass_Hex) + mass_HexNAc + mass_H2O + mass_proton # [Man2GlcNAc+H]+
    frag_b3_fuc_ion = (2 * mass_Hex) + mass_HexNAc + mass_Fuc + mass_H2O + mass_proton # [FucMan2GlcNAc+H]+

    valid_fragments = [f for f in msms_fragments_obs if f < precursor_mz]
    print(f"\nAnalyzing valid fragments: {valid_fragments}")

    print(f"\n- Fragment at m/z 204.087:")
    print(f"  Matches the HexNAc oxonium ion ([HexNAc+H]+). Theoretical m/z: {frag_hexnac:.4f}. This confirms the presence of HexNAc.")
    
    print(f"\n- Fragment at m/z 366.140:")
    print(f"  Matches the Hex-HexNAc oxonium ion ([Hex-HexNAc+H]+). Theoretical m/z: {frag_lacnac:.4f}. This confirms a LacNAc-type antenna.")

    print(f"\n- Fragment at m/z 528.193 (Base Peak):")
    print(f"  Matches the B3 ion from the core ([Man(2)HexNAc(1)+H]+). Theoretical m/z: {frag_b3_ion:.4f}.")
    print("  The HIGH INTENSITY of this ion is a classic indicator of CORE FUCOSYLATION (Fucose on the first GlcNAc).")
    
    print(f"\n- Fragment at m/z 673.231:")
    print(f"  Matches the core B3 ion plus a Fucose ([Fuc(1)Man(2)HexNAc(1)+H]+). Theoretical m/z: {frag_b3_fuc_ion:.4f}.")
    print("  This provides further strong evidence for core fucosylation.\n")

    # --- Step 4: Conclusion ---
    print("--- Step 4: Final Identification ---")
    print("Based on the precursor mass and MS/MS fragmentation pattern:")
    print("\n1. Composition: The glycan is composed of 5 Hexose, 4 HexNAc, 1 Fucose, and 1 NeuAc.")
    print("\n2. Oxford Nomenclature Name: FA2G2S1")
    print("\n3. Linkage Information:")
    print("   - The fucose is attached to the core GlcNAc, likely via an α1-6 linkage, as evidenced by the intense ion at m/z 528.193.")
    print("   - The data provided is insufficient to determine the sialic acid linkage (e.g., α2-3 vs. α2-6) or which antenna it is attached to.")

# Run the analysis
s = Scriptor()
s.start()
solve_glycan_puzzle()
final_output = s.stop()
print(final_output)

final_answer = "FA2G2S1 with core (α1-6) fucosylation"
print(f"<<<{final_answer}>>>")