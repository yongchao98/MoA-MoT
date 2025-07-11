import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# Step 1: Define all constants from the problem statement.
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol

w_oleate_large = 0.46  # mass fraction for the largest QD
w_oleate_small = 0.52  # mass fraction for the smallest QD

dH_diss_large = 70  # kJ/mol of InP for the largest QD
dH_diss_small = 120 # kJ/mol of InP for the smallest QD

dH_protonation_oleate = 7  # kJ/mol of oleate

print("--- Analysis of Enthalpy Contributions ---")
print(f"The molecular weight of InP is {MW_InP} g/mol.")
print(f"The molecular weight of oleate is {MW_oleate} g/mol.")
print(f"The enthalpy of protonation for oleate is {dH_protonation_oleate} kJ/mol.")
print("\nStep 1: Calculate the mole ratio of oleate to InP for each quantum dot size.")

# Step 2: Calculate the moles of oleate per mole of InP for each QD size.
# Mole ratio = (moles_oleate) / (moles_InP)
# Mole ratio = (mass_oleate / MW_oleate) / (mass_InP / MW_InP)
# Since mass_oleate = w_oleate and mass_InP = (1 - w_oleate) for a sample of 1g,
# Mole ratio = (w_oleate / MW_oleate) / ((1 - w_oleate) / MW_InP)

# For the largest QD
mole_ratio_large = (w_oleate_large / MW_oleate) / ((1 - w_oleate_large) / MW_InP)
print(f"Largest QD mole ratio (oleate/InP): ({w_oleate_large} / {MW_oleate}) / ((1-{w_oleate_large}) / {MW_InP}) = {mole_ratio_large:.3f}")

# For the smallest QD
mole_ratio_small = (w_oleate_small / MW_oleate) / ((1 - w_oleate_small) / MW_InP)
print(f"Smallest QD mole ratio (oleate/InP): ({w_oleate_small} / {MW_oleate}) / ((1-{w_oleate_small}) / {MW_InP}) = {mole_ratio_small:.3f}")

print("\nStep 2: Calculate the enthalpy contribution from oleate protonation (per mole of InP).")

# Step 3: Calculate the enthalpy contribution from oleate protonation for each size.
dH_contrib_large = mole_ratio_large * dH_protonation_oleate
print(f"Largest QD protonation contribution: {mole_ratio_large:.3f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_contrib_large:.2f} kJ/mol_InP")

dH_contrib_small = mole_ratio_small * dH_protonation_oleate
print(f"Smallest QD protonation contribution: {mole_ratio_small:.3f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_contrib_small:.2f} kJ/mol_InP")

print("\nStep 3: Compare the change in oleate protonation enthalpy to the total observed change.")

# Step 4: Calculate the change in oleate contribution and the total observed change.
change_in_oleate_contrib = dH_contrib_small - dH_contrib_large
print(f"Change in enthalpy due to oleate protonation = {dH_contrib_small:.2f} - {dH_contrib_large:.2f} = {change_in_oleate_contrib:.2f} kJ/mol_InP")

total_observed_change = dH_diss_small - dH_diss_large
print(f"Total observed change in dissolution enthalpy = {dH_diss_small} - {dH_diss_large} = {total_observed_change} kJ/mol_InP")

print("\n--- Conclusion ---")
print(f"The increased amount of oleate in smaller quantum dots makes the dissolution more endothermic by only {change_in_oleate_contrib:.2f} kJ/mol.")
print(f"This value is only a small fraction of the total observed increase of {total_observed_change} kJ/mol.")
print("Therefore, the protonation of oleate (Option A) cannot be the primary explanation for the large observed change in enthalpy.")
print("A much larger endothermic process that scales with the surface area must be responsible. Disrupting the tightly packed ligand shell (Option D) is the most plausible explanation for this large energy cost.")

# Restore original stdout and print captured output
sys.stdout = original_stdout
# We use a placeholder here for the final answer extraction as the prompt instructs.
# The textual explanation will guide the user to the correct choice.
output_text = captured_output.getvalue() + "\n<<<D>>>"
print(output_text)