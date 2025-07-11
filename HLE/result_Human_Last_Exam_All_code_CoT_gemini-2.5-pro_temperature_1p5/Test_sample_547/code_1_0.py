import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Constants from the problem ---
# Enthalpies in kJ/mol
dh_diss_bulk_inp = -86
dh_prot_oleate = 7
dh_diss_large_qd = 70
dh_diss_small_qd = 120

# Mass fractions (dimensionless)
mass_frac_oleate_large = 0.46
mass_frac_oleate_small = 0.52

# Molecular weights in g/mol
mw_inp = 146
mw_oleate = 281

# --- Calculations ---

# Step 1: Analyze the contribution of oleate protonation (to evaluate Choice A)

# For the largest QD
mass_frac_inp_large = 1 - mass_frac_oleate_large
# Calculate mole ratio of oleate to InP
mole_ratio_large = (mass_frac_oleate_large / mw_oleate) / (mass_frac_inp_large / mw_inp)
# Calculate the enthalpy contribution from oleate protonation per mole of InP
dh_oleate_contrib_large = mole_ratio_large * dh_prot_oleate

# For the smallest QD
mass_frac_inp_small = 1 - mass_frac_oleate_small
# Calculate mole ratio of oleate to InP
mole_ratio_small = (mass_frac_oleate_small / mw_oleate) / (mass_frac_inp_small / mw_inp)
# Calculate the enthalpy contribution from oleate protonation per mole of InP
dh_oleate_contrib_small = mole_ratio_small * dh_prot_oleate

# Step 2: Calculate the change in enthalpy contributions

# The total observed change in dissolution enthalpy
delta_dh_measured = dh_diss_small_qd - dh_diss_large_qd

# The change in enthalpy due to oleate protonation
delta_dh_oleate = dh_oleate_contrib_small - dh_oleate_contrib_large

# Step 3: Calculate the "unexplained" enthalpy (to evaluate Choice D)

# This term represents other surface effects, like ligand shell disruption
# ΔH_measured = ΔH_core + ΔH_protonation + ΔH_unexplained
# ΔH_unexplained = ΔH_measured - ΔH_core - ΔH_protonation

unknown_dh_large = dh_diss_large_qd - dh_diss_bulk_inp - dh_oleate_contrib_large
unknown_dh_small = dh_diss_small_qd - dh_diss_bulk_inp - dh_oleate_contrib_small


# --- Print the analysis ---
print("Step 1: Analyzing the effect of oleate protonation (Choice A).")
print("\nFor the largest quantum dot:")
print(f"  - The mole ratio of oleate to InP is {mole_ratio_large:.3f}.")
print(f"  - The enthalpy contribution from oleate protonation is:")
print(f"    {mole_ratio_large:.3f} (mol oleate / mol InP) * {dh_prot_oleate} (kJ / mol oleate) = {dh_oleate_contrib_large:.2f} kJ/mol of InP")

print("\nFor the smallest quantum dot:")
print(f"  - The mole ratio of oleate to InP is {mole_ratio_small:.3f}.")
print(f"  - The enthalpy contribution from oleate protonation is:")
print(f"    {mole_ratio_small:.3f} (mol oleate / mol InP) * {dh_prot_oleate} (kJ / mol oleate) = {dh_oleate_contrib_small:.2f} kJ/mol of InP")

print("\nStep 2: Comparing the change in enthalpy.")
print(f"The observed change in total enthalpy is {dh_diss_small_qd} kJ/mol - {dh_diss_large_qd} kJ/mol = {delta_dh_measured:.2f} kJ/mol.")
print(f"The change in enthalpy from just oleate protonation is {dh_oleate_contrib_small:.2f} kJ/mol - {dh_oleate_contrib_large:.2f} kJ/mol = {delta_dh_oleate:.2f} kJ/mol.")
print("Conclusion for A: This effect is far too small to explain the total observed change.")

print("\nStep 3: Calculating the large unexplained endothermic term (relevance to Choice D).")
print("This term must account for other surface effects, such as disrupting the ligand shell.")

print("\nFor the largest quantum dot, the unexplained enthalpy is:")
print(f"  ΔH_unexplained = ΔH_measured - ΔH_core - ΔH_protonation")
print(f"  ΔH_unexplained = {dh_diss_large_qd} - ({dh_diss_bulk_inp}) - {dh_oleate_contrib_large:.2f} = {unknown_dh_large:.2f} kJ/mol")

print("\nFor the smallest quantum dot, the unexplained enthalpy is:")
print(f"  ΔH_unexplained = {dh_diss_small_qd} - ({dh_diss_bulk_inp}) - {dh_oleate_contrib_small:.2f} = {unknown_dh_small:.2f} kJ/mol")

print("\nFinal Conclusion:")
print("There is a very large, positive (endothermic) enthalpy term that increases significantly with decreasing QD size.")
print("This is consistent with Choice D: the energy required to disrupt the interactions within the tightly packed ligand shell. Smaller QDs have a much larger proportion of these ligands per mole of InP, making this disruptive process more energetically costly.")

# --- Final Answer ---
# Get the content of the captured output
output_str = captured_output.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(output_str)

print("<<<D>>>")
