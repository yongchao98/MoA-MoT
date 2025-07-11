import sys
# Redirect stdout to a variable to prevent printing from the code block itself
# This is a workaround for the platform, the user should be able to run this code directly.
original_stdout = sys.stdout
from io import StringIO
sys.stdout = captured_output = StringIO()

def reaction_analysis():
    """
    Calculates the expected products from the reaction based on the stoichiometry provided
    and explains the origin of the two observed Boron NMR signals.
    """
    # --- Reaction Parameters ---
    # Molecular weight of 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    sm_mw = 317.34  # g/mol
    # Assume a starting scale for the calculation
    sm_mass = 10.0  # g

    # Equivalents of reagents used
    n_buli_eq = 1.05
    trimethyl_borate_eq = 5.0

    # --- Stoichiometric Calculations ---
    # Moles of starting material (SM)
    sm_moles = sm_mass / sm_mw

    # Moles of n-BuLi required for the desired reaction (1.0 eq)
    moles_buli_for_product = sm_moles * 1.0

    # Moles of excess n-BuLi that will form the side product (0.05 eq)
    moles_buli_for_side_product = sm_moles * (n_buli_eq - 1.0)

    # --- Output Generation ---
    print("### Analysis of Boronic Acid Synthesis ###\n")
    print(f"This script analyzes the formation of two boron-containing products in the given reaction.")
    print("-" * 50)
    print("Initial Conditions:")
    print(f"Starting Material (SM): {sm_mass:.2f} g ({sm_moles:.4f} moles) of 2-bromo-4-chloro-1-iodobenzene.")
    print(f"Reagents: {n_buli_eq} eq of n-BuLi, {trimethyl_borate_eq} eq of trimethyl borate.")
    print("-" * 50)

    print("Reaction Pathway Analysis:")
    # This represents the final desired product's mole count in the equation.
    # The numbers in the equation are: moles_buli_for_product and sm_moles
    print("1. Desired Reaction:")
    print(f"   {moles_buli_for_product:.4f} moles (1.00 eq) of n-BuLi react with {sm_moles:.4f} moles of SM.")
    print(f"   => Forms {moles_buli_for_product:.4f} moles of the desired product, (2-bromo-4-chlorophenyl)boronic acid.")
    print("\n")

    # This represents the side product's mole count in the equation.
    # The number in the equation is: moles_buli_for_side_product
    print("2. Side Reaction (The Problem):")
    print(f"   The excess {moles_buli_for_side_product:.4f} moles ({n_buli_eq - 1.0:.2f} eq) of n-BuLi remain.")
    print("   This excess n-BuLi reacts with trimethyl borate.")
    print(f"   => Forms {moles_buli_for_side_product:.4f} moles of the side product, n-butylboronic acid.")
    print("-" * 50)

    print("Conclusion:")
    print("The final reaction mixture contains two distinct boron compounds:")
    print(f"  - Desired Product: {(moles_buli_for_product / (moles_buli_for_product + moles_buli_for_side_product)) * 100:.1f}%")
    print(f"  - Side Product:    {(moles_buli_for_side_product / (moles_buli_for_product + moles_buli_for_side_product)) * 100:.1f}%")
    print("\nThe presence of these two different boronic acids explains the observation of two signals in the ¹¹B NMR spectrum.")
    print("\nProposed Solution: Use a more precise amount of n-BuLi (closer to 1.00 eq) to minimize or eliminate the formation of the n-butylboronic acid side product. This corresponds to choice C.")


# Execute the analysis
reaction_analysis()

# Restore stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())
