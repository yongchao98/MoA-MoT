import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

#--- Given Constants ---
dH_diss_bulk_inp = -86  # kJ/mol
dH_prot_oleate = 7      # kJ/mol
mw_inp = 146            # g/mol
mw_oleate = 281         # g/mol

#--- Largest Quantum Dot Data ---
dH_obs_large = 70       # kJ/mol of InP
mass_frac_oleate_large = 0.46 # 46%

#--- Smallest Quantum Dot Data ---
dH_obs_small = 120      # kJ/mol of InP
mass_frac_oleate_small = 0.52 # 52%

# --- Calculations for the LARGEST Quantum Dot ---
print("Analysis for the Largest Quantum Dot:")
# In a 100g sample
mass_oleate_large = 100 * mass_frac_oleate_large
mass_inp_large = 100 * (1 - mass_frac_oleate_large)
moles_oleate_large = mass_oleate_large / mw_oleate
moles_inp_large = mass_inp_large / mw_inp
molar_ratio_large = moles_oleate_large / moles_inp_large
print(f"Mass fraction of oleate: {mass_frac_oleate_large*100:.1f}%")
print(f"Mass fraction of InP: {(1-mass_frac_oleate_large)*100:.1f}%")
print(f"Molar ratio (mol oleate / mol InP) = ({mass_oleate_large:.1f} g / {mw_oleate} g/mol) / ({mass_inp_large:.1f} g / {mw_inp} g/mol) = {molar_ratio_large:.3f} mol/mol")

dH_oleate_contrib_large = molar_ratio_large * dH_prot_oleate
print(f"Enthalpy from oleate protonation = {molar_ratio_large:.3f} mol/mol * {dH_prot_oleate} kJ/mol = {dH_oleate_contrib_large:.2f} kJ/mol (per mol InP)")

# Total enthalpy = dH(InP diss) + dH(oleate prot) + dH(surface)
# dH(surface) = dH(observed) - dH(InP diss) - dH(oleate prot)
dH_surface_large = dH_obs_large - dH_diss_bulk_inp - dH_oleate_contrib_large
print(f"Total surface-related enthalpy = {dH_obs_large} kJ/mol - ({dH_diss_bulk_inp} kJ/mol) - {dH_oleate_contrib_large:.2f} kJ/mol = {dH_surface_large:.2f} kJ/mol")
print("-" * 30)

# --- Calculations for the SMALLEST Quantum Dot ---
print("Analysis for the Smallest Quantum Dot:")
# In a 100g sample
mass_oleate_small = 100 * mass_frac_oleate_small
mass_inp_small = 100 * (1 - mass_frac_oleate_small)
moles_oleate_small = mass_oleate_small / mw_oleate
moles_inp_small = mass_inp_small / mw_inp
molar_ratio_small = moles_oleate_small / moles_inp_small
print(f"Mass fraction of oleate: {mass_frac_oleate_small*100:.1f}%")
print(f"Mass fraction of InP: {(1-mass_frac_oleate_small)*100:.1f}%")
print(f"Molar ratio (mol oleate / mol InP) = ({mass_oleate_small:.1f} g / {mw_oleate} g/mol) / ({mass_inp_small:.1f} g / {mw_inp} g/mol) = {molar_ratio_small:.3f} mol/mol")


dH_oleate_contrib_small = molar_ratio_small * dH_prot_oleate
print(f"Enthalpy from oleate protonation = {molar_ratio_small:.3f} mol/mol * {dH_prot_oleate} kJ/mol = {dH_oleate_contrib_small:.2f} kJ/mol (per mol InP)")

dH_surface_small = dH_obs_small - dH_diss_bulk_inp - dH_oleate_contrib_small
print(f"Total surface-related enthalpy = {dH_obs_small} kJ/mol - ({dH_diss_bulk_inp} kJ/mol) - {dH_oleate_contrib_small:.2f} kJ/mol = {dH_surface_small:.2f} kJ/mol")
print("-" * 30)

# --- Conclusion ---
print("Conclusion:")
change_in_observed_dH = dH_obs_small - dH_obs_large
change_in_oleate_contrib = dH_oleate_contrib_small - dH_oleate_contrib_large
print(f"The change in enthalpy from oleate protonation ({dH_oleate_contrib_small:.2f} - {dH_oleate_contrib_large:.2f} = {change_in_oleate_contrib:.2f} kJ/mol) is far too small to explain the observed change in dissolution enthalpy ({dH_obs_small} - {dH_obs_large} = {change_in_observed_dH} kJ/mol). This rules out option A.")
print(f"The dissolution is dominated by a large endothermic term related to the surface structure, which increases from {dH_surface_large:.2f} kJ/mol to {dH_surface_small:.2f} kJ/mol as the dots get smaller.")
print("This large energy cost corresponds to overcoming the stability of the nanoparticle structure, such as breaking the bonds between the InP surface and the ligands, as well as disrupting the highly ordered, packed ligand shell. Option D provides a logical physical explanation for this large, size-dependent endothermic effect.")

# --- Final Output ---
# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output = buffer.getvalue()
# Print the content
print(output)