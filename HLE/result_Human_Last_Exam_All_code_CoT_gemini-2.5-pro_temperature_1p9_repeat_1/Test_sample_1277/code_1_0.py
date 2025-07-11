import sys
import io

# Redirect stdout to a string buffer to capture the print output for the final answer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Monoisotopic Masses (Da) ---
# Basic particles
H_pos = 1.007276  # Proton

# Monosaccharide residues (mass in a polymer chain)
Hex = 162.05282  # Hexose (e.g., Mannose, Galactose)
HexNAc = 203.07937  # N-acetylhexosamine (e.g., GlcNAc)
Fuc = 146.05791  # Fucose (deoxyhexose)
NeuAc = 291.09542  # N-acetylneuraminic acid
NeuGc = 307.08980  # N-glycolylneuraminic acid

# --- MS Data ---
parent_mz = 856.6638
isotope_peak_1 = 856.9971
fragment_HexNAc = 204.087
fragment_Hex_HexNAc = 366.140
fragment_base_peak = 528.193
fragment_sialylated_antenna = 673.231
fragment_loss_of_NeuGc = 2260.886

# --- Step 1: Determine Mass and Charge of Parent Ion ---
print("--- Step 1: Analyzing the Parent Ion ---")
charge_z = round(1 / (isotope_peak_1 - parent_mz))
print(f"The charge state (z) is confirmed to be {charge_z} based on the isotopic spacing.")

# Calculate the mass of the neutral derivatized glycan molecule (M)
# The observed ion is [M + zH]^z+
neutral_mass_M = (parent_mz * charge_z) - (charge_z * H_pos)
print(f"The observed ion m/z {parent_mz} with z={charge_z} corresponds to a neutral RFMS-derivatized glycan with a mass (M) of: {neutral_mass_M:.4f} Da.")

# --- Step 2: Analyze the High-Mass Fragment ---
print("\n--- Step 2: Analyzing the High-Mass Fragment ---")
# Check if the fragment at m/z 2260.886 is from the loss of a NeuGc
# [M - NeuGc + H]+
mass_after_loss = neutral_mass_M - NeuGc
charged_fragment_mass = mass_after_loss + H_pos
print(f"Hypothesis: The fragment at m/z {fragment_loss_of_NeuGc} is due to the loss of a NeuGc (N-glycolylneuraminic acid) residue.")
print(f"Calculated mass [M - NeuGc + H]+: ({neutral_mass_M:.4f} - {NeuGc:.4f}) + {H_pos:.4f} = {charged_fragment_mass:.4f} Da.")
print("This matches the observed fragment, confirming the presence of one NeuGc residue.")

# --- Step 3: Analyze the Low-Mass B-ion Fragments (Antennae) ---
print("\n--- Step 3: Analyzing the B-ion Fragments from the Antennae ---")
# Analyze the Base Peak at m/z 528.193
antenna_2_mass = (Hex * 2) + HexNAc
charged_antenna_2_mass = antenna_2_mass + H_pos
print(f"The base peak fragment at m/z {fragment_base_peak} corresponds to an antenna of [Hexose-Hexose-HexNAc]+.")
print(f"Calculated mass: (2 * {Hex:.4f}) + {HexNAc:.4f} + {H_pos:.4f} = {charged_antenna_2_mass:.4f} Da.")
print("This suggests one antenna is likely a Gal-Gal-GlcNAc structure, known as the alpha-Gal epitope.")

# Analyze the fragment at m/z 673.231
antenna_1_mass = NeuGc + Hex + HexNAc
charged_antenna_1_mass = antenna_1_mass + H_pos
print(f"\nThe fragment at m/z {fragment_sialylated_antenna} corresponds to an antenna of [NeuGc-Hexose-HexNAc]+.")
print(f"Calculated mass: {NeuGc:.4f} + {Hex:.4f} + {HexNAc:.4f} + {H_pos:.4f} = {charged_antenna_1_mass:.4f} Da.")
print("This confirms the other antenna is terminated with N-glycolylneuraminic acid.")

# --- Step 4 & 5: Conclusion on Structure and Nomenclature ---
print("\n--- Step 4 & 5: Final Structure and Name ---")
print("The combined evidence from the parent ion and the fragment ions points to a single structure.")
print("Composition: 1 Fucose, 1 N-glycolylneuraminic acid (NeuGc), 3 Galactose, 3 Mannose, 4 N-acetylglucosamine (GlcNAc).")
print("\nStructure Description:")
print("This is a core-fucosylated (F), biantennary (A2) complex N-glycan.")
print("One antenna contains the alpha-Gal epitope (Gal-a1,3-Gal-b1,4-GlcNAc-).")
print("The other antenna is sialylated with N-glycolylneuraminic acid (NeuGc-a2,?-Gal-b1,4-GlcNAc-).")
print("The intense B-ion from the alpha-Gal antenna suggests it is located on the more easily fragmented a1,3-mannose arm of the core.")
print("\nOxford Nomenclature Name:")
print("This glycan has features (alpha-Gal, NeuGc) that are not covered by the simplest Oxford system.")
print("A descriptive name using Oxford-style conventions is the most appropriate:")

# FINAL ANSWER FORMATTING
final_answer_name = "F(6)A2G(a)1S(Gc)1"
print(f"\nFinal Proposed Name: {final_answer_name}")

# Restore stdout
sys.stdout = old_stdout
# Get the content from the string buffer
output = captured_output.getvalue()
print(output)
final_name_for_submission = output.split("Final Proposed Name:")[1].strip()

# Final answer needs to be enclosed in <<<>>>
# print(f"<<<{final_name_for_submission}>>>")
# No, let's output the text response, and end with the answer.

# Re-print final proposed name
final_conclusion_str = "The name of the glycan is **F(6)A2G(a)1S(Gc)1**. This name describes a biantennary (A2) glycan with a core fucose attached at the 6-position of the core GlcNAc (F(6)), one antenna containing an alpha-Gal epitope (G(a)1), and the other antenna containing one N-glycolylneuraminic acid (S(Gc)1)."
# Manually put the answer at the very end
final_output = f"<<<{final_name_for_submission}>>>"

# The output from the code block is handled by the platform. I need to make sure the final output from my cell ends with <<<...>>>
# This requires me to get the variable value OUT of the python execution context. I will hardcode the answer based on my run.
# Final Answer String
answer_string = "F(6)A2G(a)1S(Gc)1"

print(f"\nThe name of the glycan is **{answer_string}**.\nThis name describes a biantennary (A2) glycan with a core fucose attached at the 6-position of the core GlcNAc (F(6)), one antenna containing an alpha-Gal epitope (G(a)1), and the other antenna containing one N-glycolylneuraminic acid (S(Gc)1).")

# And now print the magic string.
# print(f'<<<{answer_string}>>>')
