import math

def solve():
    """
    Calculates the enthalpy contribution from oleate protonation to determine
    if it explains the observed trend in dissolution enthalpy of InP quantum dots.
    """
    # Given constants
    MW_InP = 146  # g/mol
    MW_oleate = 281  # g/mol
    H_prot_oleate = 7  # kJ/mol of oleate
    H_diss_large_obs = 70  # kJ/mol of InP
    H_diss_small_obs = 120 # kJ/mol of InP

    # Data for the largest quantum dot
    mass_fraction_oleate_large = 0.46

    # Data for the smallest quantum dot
    mass_fraction_oleate_small = 0.52

    # --- Calculations for the Largest Quantum Dot ---
    # Based on a 100g sample for easier calculation
    mass_InP_large = 100.0 * (1.0 - mass_fraction_oleate_large)
    mass_oleate_large = 100.0 * mass_fraction_oleate_large

    moles_InP_large = mass_InP_large / MW_InP
    moles_oleate_large = mass_oleate_large / MW_oleate

    # Molar ratio of oleate to InP
    molar_ratio_large = moles_oleate_large / moles_InP_large
    
    # Enthalpy from oleate protonation per mole of InP
    H_prot_contribution_large = molar_ratio_large * H_prot_oleate

    # --- Calculations for the Smallest Quantum Dot ---
    # Based on a 100g sample
    mass_InP_small = 100.0 * (1.0 - mass_fraction_oleate_small)
    mass_oleate_small = 100.0 * mass_fraction_oleate_small
    
    moles_InP_small = mass_InP_small / MW_InP
    moles_oleate_small = mass_oleate_small / MW_oleate

    # Molar ratio of oleate to InP
    molar_ratio_small = moles_oleate_small / moles_InP_small

    # Enthalpy from oleate protonation per mole of InP
    H_prot_contribution_small = molar_ratio_small * H_prot_oleate

    # --- Comparison and Final Analysis ---
    observed_enthalpy_change = H_diss_small_obs - H_diss_large_obs
    explained_enthalpy_change_by_protonation = H_prot_contribution_small - H_prot_contribution_large

    print("Analysis of Enthalpy Contributions:")
    print("-" * 35)

    print("Observed Change:")
    print(f"The total dissolution enthalpy increases by {H_diss_small_obs} - {H_diss_large_obs} = {observed_enthalpy_change} kJ/mol of InP.")
    
    print("\nCalculated Contribution from Oleate Protonation:")
    # Using the equation: ΔH_change = (moles_oleate_small/moles_InP_small - moles_oleate_large/moles_InP_large) * ΔH_protonation
    print(f"For the smallest dot, the contribution is ({mass_oleate_small / MW_oleate:.3f} / {mass_InP_small / MW_InP:.3f}) * {H_prot_oleate} = {H_prot_contribution_small:.2f} kJ/mol of InP.")
    print(f"For the largest dot, the contribution is ({mass_oleate_large / MW_oleate:.3f} / {mass_InP_large / MW_InP:.3f}) * {H_prot_oleate} = {H_prot_contribution_large:.2f} kJ/mol of InP.")
    print(f"The change in enthalpy due to protonation is {H_prot_contribution_small:.2f} - {H_prot_contribution_large:.2f} = {explained_enthalpy_change_by_protonation:.2f} kJ/mol of InP.")

    print("\nConclusion:")
    print(f"The protonation of the extra oleate on the smaller dots only accounts for {explained_enthalpy_change_by_protonation:.2f} kJ/mol of the observed {observed_enthalpy_change} kJ/mol change. This is quantitatively insufficient.")
    print("Therefore, another more significant endothermic process that scales with the amount of ligand must be responsible. Option D, which attributes the endothermic trend to the energy required to disrupt the tightly packed ligand shell, is the most logical explanation. This includes breaking strong surface-ligand bonds and weaker inter-ligand interactions, both of which are more prevalent on smaller dots with higher surface area and ligand density.")

solve()