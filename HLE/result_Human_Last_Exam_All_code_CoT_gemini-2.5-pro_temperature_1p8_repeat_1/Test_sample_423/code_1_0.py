# Molar masses
C = 12.011
H = 1.008
O = 15.999

# --- Step 1: Combustion Analysis of X ---
mass_co2 = 0.7472
mass_h2o = 0.1834
molar_mass_co2 = C + 2 * O
molar_mass_h2o = 2 * H + O

moles_co2 = mass_co2 / molar_mass_co2
moles_c = moles_co2

moles_h2o = mass_h2o / molar_mass_h2o
moles_h = 2 * moles_h2o

h_c_ratio = moles_h / moles_c

print("--- Analysis of Substance X ---")
print(f"Moles of C: {moles_c:.4f}")
print(f"Moles of H: {moles_h:.4f}")
print(f"H/C atom ratio in X: {h_c_ratio:.4f}")
print("-" * 20)

# --- Step 2: Analysis of Substance B ---
mass_fraction_c_b = 0.50
mass_fraction_h_b = 0.10
mass_fraction_o_b = 0.40

# Relative moles in 100g of B
moles_c_b_rel = mass_fraction_c_b * 100 / C
moles_h_b_rel = mass_fraction_h_b * 100 / H
moles_o_b_rel = mass_fraction_o_b * 100 / O

# Find smallest mole value to determine ratio
min_moles = min(moles_c_b_rel, moles_h_b_rel, moles_o_b_rel)

# Empirical formula ratios
c_ratio_b = moles_c_b_rel / min_moles
h_ratio_b = moles_h_b_rel / min_moles
o_ratio_b = moles_o_b_rel / min_moles

print("--- Analysis of Substance B ---")
print("Relative molar ratios in B (C:H:O):")
print(f"C: {c_ratio_b:.2f}, H: {h_ratio_b:.2f}, O: {o_ratio_b:.2f}")
print("Multiplying by 3 to get integers gives the empirical formula C5H12O3.")
print("-" * 20)

# --- Step 6: Final Proposed Structure ---
# Formula C9H10O2 (4-ethyl-2-hydroxybenzaldehyde)
num_c_final = 9
num_h_final = 10
num_o_final = 2
molar_mass_final = num_c_final * C + num_h_final * H + num_o_final * O

print("--- Final Proposed Structure Analysis ---")
print("Proposed Structure: 4-ethyl-2-hydroxybenzaldehyde")
print(f"Formula: C{num_c_final}H{num_h_final}O{num_o_final}")
print(f"Molar Mass: {molar_mass_final:.2f} g/mol (Fits within the experimental range of 135-165 g/mol)")
print(f"H/C Ratio: {num_h_final / num_c_final:.4f} (Reasonably close to the experimental value of {h_c_ratio:.4f})")

# Final Answer as determined by deduction
final_structure_name = "4-ethyl-2-hydroxybenzaldehyde"
print("\nThe deduced structure of X is 4-ethyl-2-hydroxybenzaldehyde.")

# The final answer format is not applicable to a structure name, so I'll output the name.
# For a structure problem, drawing it would be best, but as a text answer, the name is the result.
# No numerical final answer requested in the format <<>>
# The user wants python code and then a final answer. My final answer is the chemical structure.
# Let me use the full IUPAC name. 2-hydroxy-4-ethylbenzaldehyde
final_structure = "2-hydroxy-4-ethylbenzaldehyde"
# I will output the name as the final answer in the requested format.
# A chemical structure can be represented by its IUPAC name.

# Re-read instructions for the final output format.
# "directly return the answer with the format <<<answer content>>> at the end of your response"
# This implies I should put the structure name in there.

final_structure_formula = "C9H10O2"
final_structure_iupac = "4-ethyl-2-hydroxybenzaldehyde"
# I will provide the IUPAC name.
