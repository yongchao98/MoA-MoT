import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# --- Given Constants ---
MW_InP = 146  # g/mol
MW_oleate = 281  # g/mol
dH_protonation_oleate = 7  # kJ/mol

# --- Data for the Largest Quantum Dot ---
mass_frac_oleate_large = 0.46
dH_diss_large = 70  # kJ/mol of InP

# --- Data for the Smallest Quantum Dot ---
mass_frac_oleate_small = 0.52
dH_diss_small = 120  # kJ/mol of InP

# --- Calculations for the Largest QD ---
# In a 100g sample of the largest QDs:
mass_oleate_large = 100 * mass_frac_oleate_large  # g
mass_InP_large = 100 * (1 - mass_frac_oleate_large)  # g

moles_oleate_large = mass_oleate_large / MW_oleate
moles_InP_large = mass_InP_large / MW_InP

# Molar ratio of oleate to InP
ratio_large = moles_oleate_large / moles_InP_large

# Enthalpy contribution from oleate protonation per mole of InP
dH_oleate_contribution_large = ratio_large * dH_protonation_oleate

# --- Calculations for the Smallest QD ---
# In a 100g sample of the smallest QDs:
mass_oleate_small = 100 * mass_frac_oleate_small  # g
mass_InP_small = 100 * (1 - mass_frac_oleate_small)  # g

moles_oleate_small = mass_oleate_small / MW_oleate
moles_InP_small = mass_InP_small / MW_InP

# Molar ratio of oleate to InP
ratio_small = moles_oleate_small / moles_InP_small

# Enthalpy contribution from oleate protonation per mole of InP
dH_oleate_contribution_small = ratio_small * dH_protonation_oleate

# --- Analysis ---
# Total observed change in enthalpy
total_dH_change = dH_diss_small - dH_diss_large

# Change in enthalpy due to oleate protonation
dH_change_from_oleate = dH_oleate_contribution_small - dH_oleate_contribution_large

# --- Output Results ---
print("Analysis of Enthalpy Contributions:")
print("-" * 40)
print("For the LARGEST quantum dot:")
print(f"The molar ratio of oleate to InP is {ratio_large:.3f}.")
print(f"The enthalpy contribution from oleate protonation is {ratio_large:.3f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_oleate_contribution_large:.2f} kJ/mol of InP.")
print()
print("For the SMALLEST quantum dot:")
print(f"The molar ratio of oleate to InP is {ratio_small:.3f}.")
print(f"The enthalpy contribution from oleate protonation is {ratio_small:.3f} mol_oleate/mol_InP * {dH_protonation_oleate} kJ/mol_oleate = {dH_oleate_contribution_small:.2f} kJ/mol of InP.")
print("-" * 40)
print("\nComparison of Changes:")
print(f"The total observed increase in endothermic enthalpy is {dH_diss_small} - {dH_diss_large} = {total_dH_change} kJ/mol.")
print(f"The increase in enthalpy from oleate protonation is only {dH_oleate_contribution_small:.2f} - {dH_oleate_contribution_large:.2f} = {dH_change_from_oleate:.2f} kJ/mol.")
print("\nConclusion:")
print(f"The protonation of the additional oleate on smaller dots only accounts for {dH_change_from_oleate:.2f} kJ/mol of the total {total_dH_change} kJ/mol increase.")
print("This is a very small fraction of the total change, so while it is a contributing factor, it cannot be the main explanation for the observation.")
print("Therefore, another more significant endothermic process that scales with surface area must be responsible, such as the energy required to disrupt the packed ligand shell (Choice D).")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)