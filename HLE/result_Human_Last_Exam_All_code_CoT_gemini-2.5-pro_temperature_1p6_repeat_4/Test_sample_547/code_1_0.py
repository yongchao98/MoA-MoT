import sys
import io

# Store the original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Constants from the problem description ---
H_diss_InP_bulk = -86  # kJ/mol
H_prot_oleate = 7      # kJ/mol
MW_InP = 146           # g/mol
MW_oleate = 281        # g/mol

# --- Data for the largest Quantum Dot ---
mass_frac_oleate_large = 0.46
mass_frac_InP_large = 1 - mass_frac_oleate_large
H_diss_exp_large = 70  # kJ/mol

# --- Data for the smallest Quantum Dot ---
mass_frac_oleate_small = 0.52
mass_frac_InP_small = 1 - mass_frac_oleate_small
H_diss_exp_small = 120 # kJ/mol

print("Evaluating the hypothesis that the total enthalpy is the sum of InP dissolution and oleate protonation (Hypothesis A).")
print("-" * 80)
print("Step 1: Calculate the molar ratio of oleate to InP for each quantum dot size.")
print("\n--- For the Largest Quantum Dot ---")
# Assume a 100g sample for calculation
moles_InP_large = (mass_frac_InP_large * 100) / MW_InP
moles_oleate_large = (mass_frac_oleate_large * 100) / MW_oleate
molar_ratio_large = moles_oleate_large / moles_InP_large
print(f"Assuming a 100g sample: Mass InP = {mass_frac_InP_large*100:.1f}g, Mass Oleate = {mass_frac_oleate_large*100:.1f}g")
print(f"Moles InP = {mass_frac_InP_large*100:.1f}g / {MW_InP} g/mol = {moles_InP_large:.4f} mol")
print(f"Moles Oleate = {mass_frac_oleate_large*100:.1f}g / {MW_oleate} g/mol = {moles_oleate_large:.4f} mol")
print(f"Molar Ratio (Oleate/InP) = {moles_oleate_large:.4f} / {moles_InP_large:.4f} = {molar_ratio_large:.4f}")

print("\n--- For the Smallest Quantum Dot ---")
moles_InP_small = (mass_frac_InP_small * 100) / MW_InP
moles_oleate_small = (mass_frac_oleate_small * 100) / MW_oleate
molar_ratio_small = moles_oleate_small / moles_InP_small
print(f"Assuming a 100g sample: Mass InP = {mass_frac_InP_small*100:.1f}g, Mass Oleate = {mass_frac_oleate_small*100:.1f}g")
print(f"Moles InP = {mass_frac_InP_small*100:.1f}g / {MW_InP} g/mol = {moles_InP_small:.4f} mol")
print(f"Moles Oleate = {mass_frac_oleate_small*100:.1f}g / {MW_oleate} g/mol = {moles_oleate_small:.4f} mol")
print(f"Molar Ratio (Oleate/InP) = {moles_oleate_small:.4f} / {moles_InP_small:.4f} = {molar_ratio_small:.4f}")

print("\n" + "-" * 80)
print("Step 2: Calculate the theoretical enthalpy of dissolution based on this model.")
# The final equation is: ΔH_calc = ΔH_InP + (Ratio * ΔH_oleate)

print("\n--- For the Largest Quantum Dot ---")
H_calc_large = H_diss_InP_bulk + (molar_ratio_large * H_prot_oleate)
print(f"Calculated ΔH = (ΔH of InP) + (Molar Ratio * ΔH of Oleate Protonation)")
print(f"Calculated ΔH = ({H_diss_InP_bulk}) + ({molar_ratio_large:.4f} * {H_prot_oleate})")
print(f"Calculated ΔH = {H_calc_large:.2f} kJ/mol")

print("\n--- For the Smallest Quantum Dot ---")
H_calc_small = H_diss_InP_bulk + (molar_ratio_small * H_prot_oleate)
print(f"Calculated ΔH = (ΔH of InP) + (Molar Ratio * ΔH of Oleate Protonation)")
print(f"Calculated ΔH = ({H_diss_InP_bulk}) + ({molar_ratio_small:.4f} * {H_prot_oleate})")
print(f"Calculated ΔH = {H_calc_small:.2f} kJ/mol")

print("\n" + "-" * 80)
print("Step 3: Compare calculated values with experimental data.")
print("\n--- For the Largest Quantum Dot ---")
print(f"Experimental ΔH = {H_diss_exp_large} kJ/mol")
print(f"Calculated ΔH (Hypothesis A) = {H_calc_large:.2f} kJ/mol")
missing_H_large = H_diss_exp_large - H_calc_large
print(f"Difference (Missing Energy) = {H_diss_exp_large} - ({H_calc_large:.2f}) = {missing_H_large:.2f} kJ/mol")


print("\n--- For the Smallest Quantum Dot ---")
print(f"Experimental ΔH = {H_diss_exp_small} kJ/mol")
print(f"Calculated ΔH (Hypothesis A) = {H_calc_small:.2f} kJ/mol")
missing_H_small = H_diss_exp_small - H_calc_small
print(f"Difference (Missing Energy) = {H_diss_exp_small} - ({H_calc_small:.2f}) = {missing_H_small:.2f} kJ/mol")

print("\n" + "=" * 80)
print("Conclusion:")
print("The calculated enthalpies based on InP dissolution and oleate protonation alone are highly exothermic (~-82 kJ/mol).")
print("The experimental enthalpies are highly endothermic (+70 to +120 kJ/mol).")
print("This shows a massive discrepancy and proves that Hypothesis A is quantitatively incorrect.")
print(f"\nThere is a large missing endothermic energy term of ~{missing_H_large:.1f} kJ/mol for the large dots and ~{missing_H_small:.1f} kJ/mol for the small dots.")
print("This missing energy is larger for the smaller dots, which have a greater proportion of ligands.")
print("This strongly supports Hypothesis D: a large amount of energy is required to disrupt the packed ligand shell, and this effect is stronger for smaller dots.")

# Restore the original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output_content = output_buffer.getvalue()

# Print the content to the real stdout
print(output_content)
# Append the final answer
print("<<<D>>>")