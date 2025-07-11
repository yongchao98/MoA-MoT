import sys
import io

# Intercept stdout to format the output as requested
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_chemistry_problem():
    """
    Calculates the outcome of a boronic acid synthesis with an inaccurate
    n-BuLi concentration to explain the observation of two boron NMR signals.
    """
    print("### Analyzing the Boronic Acid Synthesis ###\n")
    print("The goal is to synthesize (2-bromo-4-chlorophenyl)boronic acid.")
    print("The observation of two ¹¹B NMR signals suggests two boron-containing products were formed.\n")
    print("Hypothesis: An un-titrated n-BuLi solution was used, leading to a significant excess,")
    print("which created a byproduct, n-butylboronic acid, alongside the desired product.\n")
    print("--- Simulation ---")

    # --- Reaction Parameters ---
    # Molar Mass of 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    sm_mw = 317.34  # g/mol
    # Let's assume a 10 mmol reaction scale
    moles_sm = 0.010  # moles

    # --- Reagent Information ---
    # The chemist intends to add 1.05 equivalents of n-BuLi
    nBuLi_target_eq = 1.05
    # The concentration on the bottle is 1.6 M
    nBuLi_stated_M = 1.6  # mol/L
    # But due to solvent evaporation, let's assume the actual concentration is 12.5% higher
    nBuLi_actual_M = 1.8  # mol/L

    # --- Calculations ---
    # 1. Moles of n-BuLi the chemist thinks they are adding
    moles_nBuLi_intended = moles_sm * nBuLi_target_eq

    # 2. Volume of n-BuLi solution added, based on the stated concentration
    vol_nBuLi_to_add_L = moles_nBuLi_intended / nBuLi_stated_M

    # 3. Actual moles of n-BuLi delivered with that volume
    moles_nBuLi_actual = vol_nBuLi_to_add_L * nBuLi_actual_M

    # 4. Moles of n-BuLi that react with the starting material (1:1 reaction)
    moles_nBuLi_reacted_with_sm = moles_sm
    moles_desired_product = moles_nBuLi_reacted_with_sm

    # 5. Moles of n-BuLi that are in excess and form the byproduct
    moles_nBuLi_excess = moles_nBuLi_actual - moles_nBuLi_reacted_with_sm
    moles_byproduct = moles_nBuLi_excess

    print(f"Stated n-BuLi concentration: {nBuLi_stated_M:.2f} M")
    print(f"Actual n-BuLi concentration: {nBuLi_actual_M:.2f} M\n")

    print(f"Based on the stated concentration, the volume of n-BuLi added was {vol_nBuLi_to_add_L*1000:.2f} mL.")
    print(f"This delivered {moles_nBuLi_actual:.5f} moles of n-BuLi, not the intended {moles_nBuLi_intended:.5f} moles.\n")

    print("--- Final Reaction Mixture Composition ---")
    print("The final mixture contains two different boronic acids:\n")
    # Using the "output each number in the final equation" instruction format
    print(f"Equation for the final mixture:")
    print(f"Initial Reactants => {moles_desired_product:.5f} moles of Desired Product + {moles_byproduct:.5f} moles of Byproduct")
    print(f"  - Desired Product (1st B signal): (2-bromo-4-chlorophenyl)boronic acid")
    print(f"  - Byproduct (2nd B signal):   n-butylboronic acid\n")

    # --- Conclusion ---
    percent_byproduct = (moles_byproduct / (moles_desired_product + moles_byproduct)) * 100
    print(f"This creates a mixture where the byproduct is ~{percent_byproduct:.1f}% of the total boronic species.")
    print("This mixture of two products explains the two signals in the ¹¹B NMR.")
    print("\nConclusion: The problem can be solved by using a more precise amount of n-BuLi,")
    print("which requires determining its exact concentration by titration before the reaction.")


solve_chemistry_problem()

# Get the captured output and print it to the actual console
output = captured_output.getvalue()
sys.stdout = old_stdout
print(output)