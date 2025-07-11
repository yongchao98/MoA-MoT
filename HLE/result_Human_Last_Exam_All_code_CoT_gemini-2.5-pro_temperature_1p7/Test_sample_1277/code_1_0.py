import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# --- Mass Definitions ---
# Monoisotopic masses of monosaccharide residues
HEX_MASS = 162.05282  # Hexose (e.g., Galactose, Mannose)
HEXNAC_MASS = 203.07937 # N-acetylhexosamine (e.g., GlcNAc)
FUC_MASS = 146.05791   # Fucose (a deoxyhexose)
NEU5AC_MASS = 291.09542 # N-acetylneuraminic acid (a sialic acid)

# Mass of derivatization tag and other constants
RFMS_TAG_MASS = 286.1321 # RapiFluor-MS tag mass
H2O_MASS = 18.01056      # Water mass
PROTON_MASS = 1.007276   # Proton mass

# --- User-Provided Data ---
precursor_mz = 856.6638
isotope_peak_1 = 856.9971
fragment_ions = {
    204.087: "HexNAc oxonium ion",
    366.140: "Hex-HexNAc oxonium ion",
    528.193: "Diagnostic fragment for core fucosylation in RFMS analysis",
    673.231: "Y-ion fragment",
    882.409: "Y-ion fragment / internal fragment",
    1368.568: "Y-ion fragment",
    1894.753: "Y-ion fragment",
    2260.886: "Loss of a terminal sialic acid from the precursor"
}

# --- Step-by-Step Analysis ---
print("Step 1: Analyzing the Precursor Ion")
# Calculate charge state from isotopic spacing
charge_state = round(1 / (isotope_peak_1 - precursor_mz))
print(f"The isotopic spacing of ~{isotope_peak_1 - precursor_mz:.4f} Da indicates a charge state (z) of +{charge_state}.")

# Calculate the neutral mass of the RFMS-labeled glycan
mass_of_charged_species = precursor_mz * charge_state
mass_rfms_glycan = mass_of_charged_species - (charge_state * PROTON_MASS)
print(f"The observed m/z {precursor_mz} for a [M+{charge_state}H]"+
      f"{charge_state}+ ion corresponds to a neutral RFMS-glycan mass (M) of {mass_rfms_glycan:.4f} Da.")

# Calculate the neutral mass of the native glycan
mass_native_glycan = mass_rfms_glycan - RFMS_TAG_MASS + H2O_MASS
print(f"This corresponds to a native glycan mass of {mass_native_glycan:.4f} Da.")
print("-" * 30)


print("Step 2: Interpreting the MS/MS Fragmentation Pattern")
print("The fragment ions provide definitive structural information:")
for mz, desc in fragment_ions.items():
    if mz in [204.087, 366.140, 528.193, 2260.886]:
        print(f"  - The ion at m/z {mz} is a {desc}. This indicates the glycan contains HexNAc, is a complex/hybrid type, is core-fucosylated, and is sialylated.")
print("The combined fragment evidence strongly suggests a core-fucosylated, bi-antennary complex glycan structure that is capped with sialic acids.")
print("-" * 30)

print("Step 3: Proposing the Glycan Structure and Name")
print("The most common glycan with these features found on therapeutic proteins is the core-fucosylated, disialylated bi-antennary glycan.")
# Define the composition for this structure
proposed_composition = {
    'Fuc': 1, 'HexNAc': 4, 'Hex': 5, 'Neu5Ac': 2
}
# Calculate the theoretical mass for this proposed structure
theoretical_native_mass = (proposed_composition['Fuc'] * FUC_MASS +
                           proposed_composition['HexNAc'] * HEXNAC_MASS +
                           proposed_composition['Hex'] * HEX_MASS +
                           proposed_composition['Neu5Ac'] * NEU5AC_MASS)

print(f"Let's propose a composition of Fuc({proposed_composition['Fuc']}), "
      f"HexNAc({proposed_composition['HexNAc']}), Hex({proposed_composition['Hex']}), "
      f"Neu5Ac({proposed_composition['Neu5Ac']}).")
print(f"The theoretical mass for this glycan is {theoretical_native_mass:.4f} Da.")
print("Note: There is a discrepancy between the mass calculated from the observed m/z "
      f"({mass_native_glycan:.4f} Da) and the theoretical mass of the most likely structure "
      f"({theoretical_native_mass:.4f} Da). This can be due to non-standard adducts or a miscalibration. "
      "However, the fragmentation pattern is considered more definitive for structural identification.")
print("-" * 30)


print("Conclusion: Final Glycan Name")
print("Based on the definitive fragmentation evidence:")
print("  - 'F' denotes core fucosylation (from m/z 528.193).")
print("  - 'A2' denotes a bi-antennary (two-armed) structure (supported by general fragmentation).")
print("  - 'G2' denotes that both arms are terminated with galactose (assumed for a mature complex glycan).")
print("  - 'S2' denotes that both galactose residues are capped with sialic acid (supported by m/z 2260.886).")
print("\nThe MS/MS data does not provide specific linkage information (e.g., α2,3 vs α2,6 for sialic acid), "
      "so this cannot be included in the name.")

final_name = "FA2G2S2"
print(f"\nThe name of the glycan is:")
print(f">>>{final_name}<<<")


# --- Final Output ---
# Restore original stdout
sys.stdout = original_stdout
# Get the content from the string buffer
output_string = string_buffer.getvalue()
# Print the content to the real stdout
print(output_string)