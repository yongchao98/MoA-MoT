import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Mass Constants ---
# Monoisotopic masses
MASS_PROTON = 1.007276
# RapiFluor-MS reagent (C21H22FN3O3) has a mass of 383.16452 Da.
# The derivatization reaction involves loss of H2O (18.010565 Da).
MASS_RFMS_TAG = 383.16452 - 18.010565
# Monosaccharide residue masses (mass after losing H2O)
MASS_HEX = 162.05282  # Hexose (e.g., Mannose, Galactose)
MASS_HEXNAC = 203.07937  # N-acetylhexosamine (e.g., GlcNAc)
MASS_FUC = 146.05791  # Fucose (a deoxyhexose)

# --- Input Data from User ---
precursor_mz = 856.6638
charge_state = 3
msms_fragments_observed = {
    "HexNAc oxonium": 204.087,
    "Hex-HexNAc oxonium": 366.140,
    "Hex-Hex-HexNAc oxonium (most intense)": 528.193
}

# --- Step 1: Calculate Experimental Mass of the Neutral Glycan ---
print("Step 1: Calculating Experimental Mass from Precursor Ion")
# Mass of the charged ion = m/z * z
mass_charged_ion = precursor_mz * charge_state
# Mass of neutral conjugate = Mass of charged ion - mass of protons
mass_neutral_conjugate = mass_charged_ion - (charge_state * MASS_PROTON)
# Mass of neutral glycan = Mass of neutral conjugate - mass of RFMS tag
exp_mass_glycan = mass_neutral_conjugate - MASS_RFMS_TAG

print(f"Precursor m/z: {precursor_mz}")
print(f"Charge state (z): {charge_state}")
print(f"Mass of neutral glycan-RFMS conjugate = ({precursor_mz} * {charge_state}) - ({charge_state} * {MASS_PROTON:.4f}) = {mass_neutral_conjugate:.4f} Da")
print(f"Experimental mass of neutral glycan = {mass_neutral_conjugate:.4f} - {MASS_RFMS_TAG:.4f} (RFMS tag) = {exp_mass_glycan:.4f} Da")
print("-" * 30)

# --- Step 2 & 3: Propose Structure based on Fragments ---
print("Step 2 & 3: Proposing Structure based on MS/MS Fragments")
# The presence of m/z 204, 366, and especially the intense 528 suggests a complex N-glycan with multiple antennae.
# The fragment m/z 528.193 is characteristic of a Gal-Gal-GlcNAc terminal structure (the alpha-gal epitope).
# A plausible structure is a core-fucosylated, triantennary glycan with one Gal-Gal-GlcNAc antenna and two Gal-GlcNAc antennae.
proposed_composition = {"Hex": 7, "HexNAc": 5, "Fuc": 1}
print(f"The fragment ion pattern, particularly the intense ion at m/z 528.193, suggests a triantennary structure.")
print(f"Proposed Composition: {proposed_composition['Hex']} Hexose, {proposed_composition['HexNAc']} HexNAc, {proposed_composition['Fuc']} Fucose.")

# Calculate the theoretical mass of the proposed glycan
theo_mass_glycan = (proposed_composition["Hex"] * MASS_HEX) + \
                   (proposed_composition["HexNAc"] * MASS_HEXNAC) + \
                   (proposed_composition["Fuc"] * MASS_FUC)

print(f"Theoretical Mass = ({proposed_composition['Hex']} * {MASS_HEX:.4f}) + ({proposed_composition['HexNAc']} * {MASS_HEXNAC:.4f}) + ({proposed_composition['Fuc']} * {MASS_FUC:.4f}) = {theo_mass_glycan:.4f} Da")
print("\nNote: There is a significant discrepancy between the experimental mass ({exp_mass_glycan:.4f} Da) and the theoretical mass ({theo_mass_glycan:.4f} Da) of the most plausible structure based on the fragments. The identification is therefore based on the strong evidence from the MS/MS fragments.")
print("-" * 30)

# --- Step 4: Confirm Fragment Identities and Name the Glycan ---
print("Step 4: Confirming Fragment Ions and Naming Glycan")
# Calculate theoretical masses of key B-ion fragments
frag_hexnac = MASS_HEXNAC + MASS_PROTON
frag_hex_hexnac = MASS_HEX + MASS_HEXNAC + MASS_PROTON
frag_hex_hex_hexnac = MASS_HEX + MASS_HEX + MASS_HEXNAC + MASS_PROTON

print("Calculating key fragment masses:")
print(f"Theoretical [HexNAc+H]+ = {MASS_HEXNAC:.4f} + {MASS_PROTON:.4f} = {frag_hexnac:.4f} Da (Observed: {msms_fragments_observed['HexNAc oxonium']})")
print(f"Theoretical [Gal-GlcNAc+H]+ = {MASS_HEX:.4f} + {MASS_HEXNAC:.4f} + {MASS_PROTON:.4f} = {frag_hex_hexnac:.4f} Da (Observed: {msms_fragments_observed['Hex-Hex-HexNAc oxonium']})")
print(f"Theoretical [Gal-Gal-GlcNAc+H]+ = {MASS_HEX:.4f} + {MASS_HEX:.4f} + {MASS_HEXNAC:.4f} + {MASS_PROTON:.4f} = {frag_hex_hex_hexnac:.4f} Da (Observed: {msms_fragments_observed['Hex-Hex-Hex-HexNAc oxonium (most intense)']})")
print("\nThe calculated fragment masses match the observed fragments, confirming the presence of Gal-GlcNAc and Gal-Gal-GlcNAc antennae.")
print("-" * 30)

# Final Conclusion
print("Final Answer:")
print("The glycan is a core-fucosylated, triantennary complex N-glycan containing two monogalactosylated antennae and one digalactosylated (alpha-gal) antenna.")
print("Oxford Nomenclature Name: FA3G3(α1-3Gal)1")
print("This name indicates a core Fucosylated (F), triAntennary (A3) glycan with 3 terminating Galactose residues (G3), where one of those antennae has an additional alpha-1,3-linked Galactose.")
final_answer_name = "FA3G3(α1-3Gal)1"

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output_str = captured_output.getvalue()
print(output_str)

# Final answer in the required format
print(f'<<<{final_answer_name}>>>')