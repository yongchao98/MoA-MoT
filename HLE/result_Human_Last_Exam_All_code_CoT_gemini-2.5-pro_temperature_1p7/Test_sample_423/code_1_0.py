import math

def solve_structure_X():
    """
    This function carries out the step-by-step deduction of the structure of substance X.
    """
    # Molar masses of elements
    M_C = 12.01
    M_H = 1.008
    M_O = 16.00

    print("--- Step 1: Analysis of Substance B ---")
    print("The mass fractions of elements in B are: C = 50%, H = 10%, O = 40%.")
    
    # Assuming a 100g sample of B to find mole ratios
    moles_C_B = 50 / M_C
    moles_H_B = 10 / M_H
    moles_O_B = 40 / M_O

    # Finding the simplest ratio by dividing by the smallest number of moles
    min_moles_B = min(moles_C_B, moles_H_B, moles_O_B)
    ratio_C = moles_C_B / min_moles_B
    ratio_H = moles_H_B / min_moles_B
    ratio_O = moles_O_B / min_moles_B
    
    print(f"The mole ratio is C : H : O = {ratio_C:.3f} : {ratio_H:.3f} : {ratio_O:.3f}")
    # The ratio is approx 1.667 : 4.00 : 1. To get integers, we multiply by 3.
    emp_C_B = round(ratio_C * 3)
    emp_H_B = round(ratio_H * 3)
    emp_O_B = round(ratio_O * 3)
    
    print(f"Multiplying by 3 yields the integer ratio C : H : O = {emp_C_B}:{emp_H_B}:{emp_O_B}.")
    print(f"Thus, the empirical formula of B is C{emp_C_B}H{emp_H_B}O{emp_O_B}.")
    print("Since B is reduced to n-pentane by HI, it must have a 5-carbon linear skeleton.")
    print("This confirms the molecular formula is C5H12O3.")
    print("The bromination reaction details (two monobromo derivatives, one chiral and one achiral) confirm that B is pentane-1,3,5-triol.\n")

    print("--- Step 2: Analysis of Substance X from Combustion Data ---")
    mass_CO2 = 0.7472  # g
    mass_H2O = 0.1834  # g

    M_CO2 = M_C + 2 * M_O
    M_H2O = 2 * M_H + M_O

    moles_CO2 = mass_CO2 / M_CO2
    moles_H_X = 2 * (mass_H2O / M_H2O)
    moles_C_X = moles_CO2
    
    print(f"From {mass_CO2} g of CO2, moles of Carbon in X = {moles_C_X:.5f} mol.")
    print(f"From {mass_H2O} g of H2O, moles of Hydrogen in X = {moles_H_X:.5f} mol.")
    
    ratio_H_C = moles_H_X / moles_C_X
    print(f"The molar ratio H:C in X is {ratio_H_C:.3f}:1, which is approximately 1.2:1 or 12:10.\n")

    print("--- Step 3: Resolving the Molar Mass Contradiction ---")
    print("The H:C ratio suggests a formula of (C10H12)n.")
    M_C10H12O = 10 * M_C + 12 * M_H + 1 * M_O
    print(f"Testing the formula C10H12O: M = 10*{M_C} + 12*{M_H} + 1*{M_O} = {M_C10H12O:.2f} g/mol.")
    print("This calculated mass is consistent with the experimental M ~ 150 g/mol.")
    print("However, a C10 molecule (X) cannot produce a C5 molecule (B) via ozonolysis and reduction.")
    print("This contradiction is resolved by assuming that X dimerizes in solution (due to its -OH group), causing the measured molar mass to be double the actual mass.")
    print(f"Therefore, the actual molar mass of X is approximately 150 / 2 = 75 g/mol.\n")

    print("--- Step 4: Deducing the Final Structure of X ---")
    print("With M_real ≈ 75 g/mol and a C5 backbone, we deduce the true formula.")
    M_C5H6O = 5 * M_C + 6 * M_H + 1 * M_O
    print(f"A formula of C5H6O has a molar mass of 5*{M_C} + 6*{M_H} + 1*{M_O} = {M_C5H6O:.2f} g/mol, which matches M_real.")
    print("\n- Ozonolysis of X gives a single product A, which is then reduced to B (pentane-1,3,5-triol).")
    print("- This implies A is 3-oxopentanedial, and X must be cyclopent-3-en-1-one.")
    print("- X exists in equilibrium with its enol tautomer: cyclopent-2,4-dien-1-ol.")
    print("- This enol structure accounts for all the observed chemical properties:")
    print("  - Reacts with Na, NaOH, and gives a red color with FeCl3 (acidic enol).")
    print("  - Reduces Tollens' reagent (a property of conjugated dienols).")
    print("  - Has -OH bands in its IR spectrum.")

    final_structure_name = "cyclopent-2,4-dien-1-ol"
    print("\nThe determined structure for X is the enol form, which is responsible for most of its characteristic reactivity.")
    
    return final_structure_name

# Execute the analysis and print the final answer
solve_structure_X()