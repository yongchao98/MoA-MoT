import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Given Constants ---
h_diss_bulk_inp = -86  # kJ/mol
h_prot_oleate = 7      # kJ/mol
h_diss_large_qd = 70   # kJ/mol of InP
h_diss_small_qd = 120  # kJ/mol of InP
mass_frac_oleate_large = 0.46
mass_frac_oleate_small = 0.52
mw_inp = 146  # g/mol
mw_oleate = 281 # g/mol

# --- Calculations for the LARGEST Quantum Dot ---
# Based on a 100g sample
mass_inp_large = 100 * (1 - mass_frac_oleate_large)
moles_oleate_large = (100 * mass_frac_oleate_large) / mw_oleate
moles_inp_large = mass_inp_large / mw_inp
ratio_large = moles_oleate_large / moles_inp_large
h_oleate_contrib_large = ratio_large * h_prot_oleate

# --- Calculations for the SMALLEST Quantum Dot ---
# Based on a 100g sample
mass_inp_small = 100 * (1 - mass_frac_oleate_small)
moles_oleate_small = (100 * mass_frac_oleate_small) / mw_oleate
moles_inp_small = mass_inp_small / mw_inp
ratio_small = moles_oleate_small / moles_inp_small
h_oleate_contrib_small = ratio_small * h_prot_oleate

# --- Analysis of Enthalpy Changes ---
delta_h_observed = h_diss_small_qd - h_diss_large_qd
delta_h_oleate = h_oleate_contrib_small - h_oleate_contrib_large

# --- Calculating the "Surface Penalty" Term ---
# H_penalty = H_measured - H_bulk_diss - H_protonation
h_penalty_large = h_diss_large_qd - h_diss_bulk_inp - h_oleate_contrib_large
h_penalty_small = h_diss_small_qd - h_diss_bulk_inp - h_oleate_contrib_small

# --- Output the analysis ---
print("--- Analysis of Energy Contributions ---")
print("\n1. Contribution from Oleate Protonation (evaluating Choice A):")
print(f"The increase in enthalpy from oleate protonation when going from large to small QDs is only {delta_h_oleate:.2f} kJ/mol.")
print(f"The total observed enthalpy increase is {delta_h_observed} kJ/mol.")
print("Conclusion: The contribution from oleate protonation is far too small to explain the total observed change.")

print("\n2. Calculation of the Full 'Surface Penalty':")
print("The overall dissolution enthalpy equation can be arranged to solve for the surface penalty:")
print("H_penalty = H_measured - H_diss(bulk) - H_protonation(oleate)")

print("\nFor the LARGEST quantum dot:")
print(f"H_penalty = {h_diss_large_qd} kJ/mol - ({h_diss_bulk_inp} kJ/mol) - {h_oleate_contrib_large:.2f} kJ/mol")
print(f"Resulting H_penalty = {h_penalty_large:.2f} kJ/mol")

print("\nFor the SMALLEST quantum dot:")
print(f"H_penalty = {h_diss_small_qd} kJ/mol - ({h_diss_bulk_inp} kJ/mol) - {h_oleate_contrib_small:.2f} kJ/mol")
print(f"Resulting H_penalty = {h_penalty_small:.2f} kJ/mol")

print("\n--- Final Conclusion ---")
print("The calculations show a large, positive 'surface penalty' that increases significantly (by nearly 50 kJ/mol) as the quantum dot size decreases.")
print("This large endothermic term, which scales with the surface-to-volume ratio, is best explained by the energy required to disrupt the tightly packed organic ligand shell on the nanoparticle surface.")
print("This matches the description in Choice D.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()
print(output_string)