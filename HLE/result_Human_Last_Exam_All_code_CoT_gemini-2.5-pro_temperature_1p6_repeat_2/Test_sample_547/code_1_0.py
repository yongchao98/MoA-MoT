import math

def analyze_quantum_dot_dissolution():
    """
    Analyzes the enthalpy changes during the dissolution of InP quantum dots
    to determine the primary cause of the observed trend.
    """
    # Given constants
    MW_InP = 146  # g/mol
    MW_oleate = 281  # g/mol
    H_protonation_oleate = 7  # kJ/mol of oleate

    # --- Data for the largest quantum dot ---
    mass_fraction_oleate_large = 0.46
    H_dissolution_large = 70  # kJ/mol of InP

    # --- Data for the smallest quantum dot ---
    mass_fraction_oleate_small = 0.52
    H_dissolution_small = 120  # kJ/mol of InP

    # --- Calculations for the largest QD (based on a 100g sample) ---
    # Mass of each component in 100g of QD material
    mass_InP_large = 100 * (1 - mass_fraction_oleate_large)
    mass_oleate_large = 100 * mass_fraction_oleate_large

    # Moles of each component
    moles_InP_large = mass_InP_large / MW_InP
    moles_oleate_large = mass_oleate_large / MW_oleate

    # Molar ratio of oleate to InP
    ratio_large = moles_oleate_large / moles_InP_large

    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_from_oleate_large = ratio_large * H_protonation_oleate

    # --- Calculations for the smallest QD (based on a 100g sample) ---
    # Mass of each component in 100g of QD material
    mass_InP_small = 100 * (1 - mass_fraction_oleate_small)
    mass_oleate_small = 100 * mass_fraction_oleate_small

    # Moles of each component
    moles_InP_small = mass_InP_small / MW_InP
    moles_oleate_small = mass_oleate_small / MW_oleate

    # Molar ratio of oleate to InP
    ratio_small = moles_oleate_small / moles_InP_small

    # Enthalpy contribution from oleate protonation (per mole of InP)
    h_from_oleate_small = ratio_small * H_protonation_oleate

    # --- Analysis of the change in enthalpy ---
    # Total observed change in dissolution enthalpy
    delta_H_observed = H_dissolution_small - H_dissolution_large

    # Change in enthalpy that can be explained by oleate protonation
    delta_H_from_oleate = h_from_oleate_small - h_from_oleate_large
    
    # The remaining change must be due to other factors (surface energy, ligand shell disruption, etc.)
    delta_H_other = delta_H_observed - delta_H_from_oleate

    print("--- Analysis of Enthalpy Contributions ---")
    print(f"For the largest QD, the enthalpy from oleate protonation is {h_from_oleate_large:.2f} kJ/mol of InP.")
    print(f"For the smallest QD, the enthalpy from oleate protonation is {h_from_oleate_small:.2f} kJ/mol of InP.")
    print("\n--- Breakdown of the Enthalpy Change ---")
    print("The total observed increase in enthalpy from largest to smallest QDs is:")
    print(f"ΔH_observed = {H_dissolution_small} - {H_dissolution_large} = {delta_H_observed:.2f} kJ/mol")
    
    print("\nThe portion of this increase due to more oleate being protonated is:")
    print(f"ΔH_oleate = {h_from_oleate_small:.2f} - {h_from_oleate_large:.2f} = {delta_H_from_oleate:.2f} kJ/mol")
    
    print("\nThis means the majority of the change is from other factors related to the QD structure:")
    print("ΔH_other = ΔH_observed - ΔH_oleate")
    print("Final Equation:")
    print(f"{delta_H_observed:.2f} = {delta_H_from_oleate:.2f} + {delta_H_other:.2f}")

    print("\n--- Conclusion ---")
    print("The calculation shows that the increased protonation of oleate only accounts for about 0.84 kJ/mol of the 50.00 kJ/mol observed increase in endothermicity. Therefore, answer choice A is incorrect.")
    print("The remaining 49.16 kJ/mol must come from another source related to the increased surface area of smaller dots. Option D provides the most logical physical explanation: disrupting the tightly packed shell of organic ligands, including breaking ligand-nanocrystal bonds and inter-ligand interactions, is a highly endothermic process that is much more significant for smaller quantum dots, which have a far greater proportion of surface ligands.")

analyze_quantum_dot_dissolution()