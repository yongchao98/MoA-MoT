import sys

def solve_enthalpy_problem():
    """
    Analyzes the enthalpy changes during the dissolution of InP quantum dots
    to determine the most logical explanation for the observed phenomena.
    """

    # --- Given Constants ---
    MW_InP = 146  # g/mol
    MW_oleate = 281  # g/mol

    dH_diss_bulk_InP = -86  # kJ/mol
    dH_protonation_oleate = 7  # kJ/mol

    # Largest Quantum Dot Data
    mass_frac_oleate_large = 0.46
    dH_diss_large_QD = 70  # kJ/mol

    # Smallest Quantum Dot Data
    mass_frac_oleate_small = 0.52
    dH_diss_small_QD = 120  # kJ/mol

    # --- Step 1: Analyze the Largest Quantum Dot ---
    print("--- Analysis for the Largest Quantum Dot ---")
    # In 100g of QD material:
    mass_oleate_large = 100 * mass_frac_oleate_large
    mass_InP_large = 100 - mass_oleate_large
    moles_oleate_large = mass_oleate_large / MW_oleate
    moles_InP_large = mass_InP_large / MW_InP
    mole_ratio_large = moles_oleate_large / moles_InP_large
    
    print(f"Molar Ratio (Oleate/InP): ({mass_oleate_large:.1f}g / {MW_oleate}g/mol) / ({mass_InP_large:.1f}g / {MW_InP}g/mol) = {mole_ratio_large:.4f} mol oleate / mol InP")

    # --- Step 2: Calculate Enthalpy Contribution from Oleate Protonation (Largest QD) ---
    dH_protonation_contrib_large = mole_ratio_large * dH_protonation_oleate
    print("Enthalpy from Oleate Protonation (per mole InP):")
    print(f"Equation: (moles_oleate / moles_InP) * ΔH_protonation_oleate")
    print(f"Calculation: {mole_ratio_large:.4f} * {dH_protonation_oleate} kJ/mol = {dH_protonation_contrib_large:.2f} kJ/mol\n")


    # --- Step 3: Analyze the Smallest Quantum Dot ---
    print("--- Analysis for the Smallest Quantum Dot ---")
    # In 100g of QD material:
    mass_oleate_small = 100 * mass_frac_oleate_small
    mass_InP_small = 100 - mass_oleate_small
    moles_oleate_small = mass_oleate_small / MW_oleate
    moles_InP_small = mass_InP_small / MW_InP
    mole_ratio_small = moles_oleate_small / moles_InP_small
    
    print(f"Molar Ratio (Oleate/InP): ({mass_oleate_small:.1f}g / {MW_oleate}g/mol) / ({mass_InP_small:.1f}g / {MW_InP}g/mol) = {mole_ratio_small:.4f} mol oleate / mol InP")
    
    # --- Step 4: Calculate Enthalpy Contribution from Oleate Protonation (Smallest QD) ---
    dH_protonation_contrib_small = mole_ratio_small * dH_protonation_oleate
    print("Enthalpy from Oleate Protonation (per mole InP):")
    print(f"Equation: (moles_oleate / moles_InP) * ΔH_protonation_oleate")
    print(f"Calculation: {mole_ratio_small:.4f} * {dH_protonation_oleate} kJ/mol = {dH_protonation_contrib_small:.2f} kJ/mol\n")

    # --- Step 5: Compare Enthalpy Changes (Test of Choice A) ---
    print("--- Comparing Enthalpy Changes ---")
    total_dH_change = dH_diss_small_QD - dH_diss_large_QD
    protonation_dH_change = dH_protonation_contrib_small - dH_protonation_contrib_large
    
    print(f"Total observed change in dissolution enthalpy = {dH_diss_small_QD} - {dH_diss_large_QD} = {total_dH_change} kJ/mol")
    print(f"Change due to oleate protonation = {dH_protonation_contrib_small:.2f} - {dH_protonation_contrib_large:.2f} = {protonation_dH_change:.2f} kJ/mol")
    print("Observation: The change due to protonation accounts for only a tiny fraction of the total change. Therefore, explanation A is not sufficient.\n")
    
    # --- Step 6: Calculate the 'Surface Enthalpy' Term ---
    # ΔH_surface = ΔH_diss(QD) - ΔH_diss(bulk) - ΔH_protonation_contribution
    print("--- Calculating the Excess Surface Enthalpy ---")
    
    # Largest QD
    dH_surface_large = dH_diss_large_QD - dH_diss_bulk_InP - dH_protonation_contrib_large
    print("Surface Enthalpy (Largest QD) = ΔH_diss(QD) - ΔH_diss(bulk) - ΔH_protonation")
    print(f"Calculation: {dH_diss_large_QD} kJ/mol - ({dH_diss_bulk_InP} kJ/mol) - {dH_protonation_contrib_large:.2f} kJ/mol = {dH_surface_large:.2f} kJ/mol")

    # Smallest QD
    dH_surface_small = dH_diss_small_QD - dH_diss_bulk_InP - dH_protonation_contrib_small
    print("Surface Enthalpy (Smallest QD) = ΔH_diss(QD) - ΔH_diss(bulk) - ΔH_protonation")
    print(f"Calculation: {dH_diss_small_QD} kJ/mol - ({dH_diss_bulk_InP} kJ/mol) - {dH_protonation_contrib_small:.2f} kJ/mol = {dH_surface_small:.2f} kJ/mol\n")

    print("Conclusion:")
    print("The calculations show that a very large, positive 'surface enthalpy' term is required to explain the observations.")
    print(f"This term increases from {dH_surface_large:.2f} kJ/mol for the largest dots to {dH_surface_small:.2f} kJ/mol for the smallest dots.")
    print("This large energy cost is associated with breaking down the highly structured surface of the quantum dot, including disrupting the tightly packed ligand shell.")
    print("As smaller dots have a greater surface-area-to-volume ratio, this surface effect becomes more pronounced, making their dissolution more endothermic.")
    print("This strongly supports choice D.")

# Run the analysis
solve_enthalpy_problem()
# Redirect stdout to a variable to prevent final answer from being printed here
original_stdout = sys.stdout
sys.stdout = None
print("<<<D>>>")
sys.stdout = original_stdout
