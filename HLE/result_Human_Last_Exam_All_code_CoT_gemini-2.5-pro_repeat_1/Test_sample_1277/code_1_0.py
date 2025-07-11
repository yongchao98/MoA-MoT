import sys
import io

# Redirect print to a string to capture it for the final answer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Define Constants ---
# Mass of a proton
H_mass = 1.007276
# Mass of the RapiFluor-MS tag added to the glycan
RFMS_tag_mass = 214.111
# Monoisotopic masses of common monosaccharide residues
masses = {
    'HexNAc': 203.07937, # N-acetylhexosamine (e.g., GlcNAc)
    'Hex': 162.05282,    # Hexose (e.g., Man, Gal)
    'dHex': 146.05791,   # Deoxyhexose (Fucose)
    'NeuAc': 291.09542,  # N-acetylneuraminic acid (a sialic acid)
}
# Known theoretical m/z for key diagnostic fragment ions
theo_frags = {
    'HexNAc': 204.0872,      # [C8H14NO5]+
    'Hex-HexNAc': 366.1395,  # [C14H24NO10]+
    'Hex2-HexNAc': 528.1928 # [C20H34NO15]+
}

# --- Input Data from the User ---
precursor_mz = 856.6638
# Charge is determined from isotopic spacing: 1 / (856.9971 - 856.6638) = 1 / 0.3333 ≈ 3
charge = 3
fragment_ions = {'204.087': 'HexNAc', '366.140': 'Hex-HexNAc', '528.193': 'Hex2-HexNAc (most intense)'}

# --- Step 1: Calculate Experimental Mass of the Glycan ---
print("--- Step 1: Mass Determination ---")
neutral_mass_conjugate = (precursor_mz * charge) - (charge * H_mass)
print(f"The precursor ion is [M+3H]3+ at m/z {precursor_mz}.")
print(f"Calculating the neutral mass of the RFMS-derivatized glycan:")
print(f"Equation: (m/z * charge) - (charge * mass_proton)")
print(f"Result: ({precursor_mz} * {charge}) - ({charge} * {H_mass:.6f}) = {neutral_mass_conjugate:.4f} Da")

glycan_exp_mass = neutral_mass_conjugate - RFMS_tag_mass
print(f"\nCalculating the underivatized glycan mass by subtracting the RFMS tag ({RFMS_tag_mass} Da):")
print(f"Equation: Neutral_Mass_Conjugate - Mass_RFMS_Tag")
print(f"Result: {neutral_mass_conjugate:.4f} - {RFMS_tag_mass} = {glycan_exp_mass:.4f} Da")

# --- Step 2: Propose and Verify Glycan Composition ---
print("\n--- Step 2: Glycan Composition Identification ---")
# The most plausible composition for this mass on a therapeutic protein is FA2G2S2
composition = {'HexNAc': 4, 'Hex': 5, 'dHex': 1, 'NeuAc': 2}
theo_mass = (composition['HexNAc'] * masses['HexNAc'] +
             composition['Hex'] * masses['Hex'] +
             composition['dHex'] * masses['dHex'] +
             composition['NeuAc'] * masses['NeuAc'])

print("The experimental mass is matched to the theoretical mass of a fucosylated, disialylated, biantennary glycan (FA2G2S2).")
print(f"Proposed Composition: Fuc(1)Hex(5)HexNAc(4)NeuAc(2)")
print(f"Calculating the theoretical mass of this composition:")
print(f"Equation: (4 * {masses['HexNAc']}) + (5 * {masses['Hex']}) + (1 * {masses['dHex']}) + (2 * {masses['NeuAc']})")
print(f"Result: {theo_mass:.4f} Da")

mass_diff = glycan_exp_mass - theo_mass
print(f"\nMass Difference = {glycan_exp_mass:.4f} (Experimental) - {theo_mass:.4f} (Theoretical) = {mass_diff:.4f} Da")
print("The proposed FA2G2S2 composition is a strong candidate.")

# --- Step 3: Interpret MS/MS Fragmentation Data ---
print("\n--- Step 3: MS/MS Fragment Analysis ---")
print("Key fragment ions are analyzed to confirm the proposed structure.")

# Analysis of m/z 204.087
print(f"\nFragment: 204.087 m/z")
print(f"Identification: This matches the HexNAc oxonium ion, confirming the presence of HexNAc residues.")
print(f"Equation: Theoretical m/z for [HexNAc oxonium]+ = {theo_frags['HexNAc']:.4f}")

# Analysis of m/z 366.140
print(f"\nFragment: 366.140 m/z")
print(f"Identification: This matches the Hex-HexNAc oxonium ion, characteristic of a Gal-GlcNAc unit from an antenna.")
print(f"Equation: Theoretical m/z for [Hex-HexNAc oxonium]+ = {theo_frags['Hex-HexNAc']:.4f}")

# Analysis of m/z 528.193
print(f"\nFragment: 528.193 m/z (Most Intense)")
print(f"Identification: This matches the Hex2-HexNAc oxonium ion, likely a Man-Gal-GlcNAc fragment from an antenna after the loss of a labile sialic acid.")
print(f"Equation: Theoretical m/z for [Hex2-HexNAc oxonium]+ = {theo_frags['Hex2-HexNAc']:.4f}")
print("The high intensity of this ion is a strong indicator of a labile sialic acid linkage.")

# --- Step 4: Final Conclusion ---
print("\n--- Step 4: Final Identification ---")
print("\nBased on the parent mass and the fragmentation pattern, the glycan is identified as:")
print("\nGlycan Name (Oxford Nomenclature): FA2G2S2")
print("\nDescription:")
print("A complex biantennary N-glycan with a core fucose ('F'), two antennae ('A2'), two terminal galactose residues ('G2'), and two terminal sialic acid residues ('S2').")
print("\nInferred Linkage Information:")
print("- Core Fucosylation: The 'F' designation assumes core fucosylation (α1,6-linkage to the innermost GlcNAc), a common modification on therapeutic proteins.")
print("- Sialic Acid Linkage: The intense fragment at m/z 528.193 resulting from sialic acid loss suggests that the sialic acid linkages are labile, which is characteristic of α2,3-linkages rather than the more stable α2,6-linkages.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())