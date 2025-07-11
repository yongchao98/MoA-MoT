import math

def solve_boronic_acid_problem():
    """
    Explains the cause of observing two Boron NMR signals and calculates the
    stoichiometry to demonstrate the solution.
    """

    # --- 1. Define Reaction Parameters ---
    # Starting Material (SM): 2-bromo-4-chloro-1-iodobenzene
    sm_mw = 317.34  # g/mol

    # Reagents
    nBuLi_eq = 1.05  # Equivalents of n-BuLi
    tmb_eq = 5.0  # Equivalents of trimethyl borate

    # --- 2. Explain the Core Problem ---
    print("--- Problem Analysis ---")
    print("The observation of two ¹¹B NMR signals suggests a mixture of products.")
    print("A likely scenario is the formation of two species:")
    print("  1. Desired Product: (2-bromo-4-chlorophenyl)boronic acid, Ar-B(OH)₂")
    print("  2. Side Product:   Diarylborinic acid, Ar₂B-OH\n")

    print("--- Cause of Side Product ---")
    print("The diarylborinic acid (Ar₂B-OH) is formed when two equivalents of the")
    print("aryllithium intermediate (Ar-Li) react with one molecule of trimethyl borate.")
    print("This happens when there is an excess of Ar-Li in the reaction mixture.\n")

    print("A primary reason for excess Ar-Li is an inaccurate concentration of the")
    print("n-BuLi solution. The concentration of commercial n-BuLi can change over")
    print("time due to degradation or solvent evaporation. Relying on the label's")
    print("concentration can lead to significant stoichiometric errors.\n")

    # --- 3. Stoichiometric Calculation Example ---
    print("--- Illustrative Calculation ---")
    scale_mmol = 20.0  # Let's assume a 20.0 mmol scale reaction
    sm_mass = (scale_mmol / 1000) * sm_mw

    print(f"Reaction Scale:")
    print(f"  - Starting Material (SM) moles = {scale_mmol:.1f} mmol")
    print(f"  - Starting Material (SM) mass = {sm_mass:.2f} g\n")

    # Moles of n-BuLi needed based on desired equivalents
    nBuLi_moles_needed = scale_mmol * nBuLi_eq
    
    print("Stoichiometric Equation for n-BuLi:")
    print(f"  Moles n-BuLi = Moles SM * Equivalents")
    print(f"  Moles n-BuLi = {scale_mmol:.1f} mmol * {nBuLi_eq}")
    print(f"  Calculated Moles n-BuLi = {nBuLi_moles_needed:.2f} mmol\n")

    # Demonstrate the error
    labeled_nBuLi_conc = 2.5  # M (mol/L)
    actual_nBuLi_conc = 2.8  # M (e.g., after some solvent evaporation, confirmed by titration)

    volume_to_add_based_on_label = (nBuLi_moles_needed / 1000) / labeled_nBuLi_conc * 1000  # in mL
    
    print("Scenario: Using Inaccurate Labeled Concentration")
    print(f"  - Labeled n-BuLi concentration = {labeled_nBuLi_conc:.2f} M")
    print(f"  - Volume added based on label = {nBuLi_moles_needed:.2f} mmol / {labeled_nBuLi_conc:.2f} M = {volume_to_add_based_on_label:.2f} mL")

    # Actual moles added if concentration is different
    actual_moles_added = (volume_to_add_based_on_label / 1000) * actual_nBuLi_conc * 1000 # in mmol
    excess_moles = actual_moles_added - scale_mmol # Excess over the starting material

    print(f"  - Actual n-BuLi concentration (from titration) = {actual_nBuLi_conc:.2f} M")
    print(f"  - Moles actually added = {volume_to_add_based_on_label:.2f} mL * {actual_nBuLi_conc:.2f} M = {actual_moles_added:.2f} mmol\n")
    
    print("Final Equation for Excess Reagent:")
    print(f"  Excess Ar-Li (mmol) = Actual Moles Added - Moles SM")
    # Here we output each number in the final equation
    print(f"  Excess Ar-Li (mmol) = {actual_moles_added:.2f} - {scale_mmol:.1f} = {excess_moles:.2f} mmol")
    print(f"This {excess_moles:.2f} mmol of excess Ar-Li will react with the boronic ester product,")
    print("leading to the diarylborinic acid impurity.\n")
    
    # --- 4. Conclusion ---
    print("--- Solution ---")
    print("To avoid the side reaction, the n-BuLi solution must be titrated immediately")
    print("before use to determine its exact concentration. Then, the precise volume")
    print("for the desired 1.05 equivalents can be added.")
    print("Therefore, the best way to solve the problem is to use a more precise amount of n-BuLi.")

if __name__ == "__main__":
    solve_boronic_acid_problem()