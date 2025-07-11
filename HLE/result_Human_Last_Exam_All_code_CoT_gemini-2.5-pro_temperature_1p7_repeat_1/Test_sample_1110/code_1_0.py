import sys

def solve_chemistry_problem():
    """
    Calculates the precise volume of n-BuLi needed for a selective monolithiation reaction.

    The key to avoiding the formation of a di-boronated side product is to use a precise
    stoichiometric amount of n-BuLi. This requires knowing the exact concentration
    of the n-BuLi solution via titration.
    """

    # Atomic masses in g/mol
    ATOMIC_MASSES = {
        'C': 12.011,
        'H': 1.008,
        'Br': 79.904,
        'Cl': 35.453,
        'I': 126.90
    }

    # --- Step 0: Define experiment parameters ---
    # These values should be adjusted for your specific experiment.
    # Mass of the starting material, 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    mass_starting_material = 5.0  # grams

    # CONCENTRATION OF n-BuLi AFTER TITRATION (this is the crucial, measured value)
    nBuLi_concentration_M = 1.55  # mol/L (M)

    # To ensure selectivity, we use exactly 1.00 equivalent.
    nBuLi_equivalents = 1.00

    print("Plan: Calculate the precise volume of n-BuLi to add for a 1:1 reaction.\n")

    # --- Step 1: Calculate the Molar Mass of the starting material ---
    # Formula: C6H3BrClI
    molar_mass_sm = (6 * ATOMIC_MASSES['C'] +
                     3 * ATOMIC_MASSES['H'] +
                     1 * ATOMIC_MASSES['Br'] +
                     1 * ATOMIC_MASSES['Cl'] +
                     1 * ATOMIC_MASSES['I'])
    print(f"Step 1: Calculated Molar Mass of C6H3BrClI = {molar_mass_sm:.2f} g/mol")
    print("-" * 50)

    # --- Step 2: Calculate moles of the starting material ---
    # Equation: Moles = Mass / Molar Mass
    moles_sm = mass_starting_material / molar_mass_sm
    print("Step 2: Calculate moles of starting material.")
    print(f"  Equation: Moles = Mass / Molar Mass")
    print(f"  Calculation: {mass_starting_material:.2f} g / {molar_mass_sm:.2f} g/mol")
    print(f"  Result: {moles_sm:.6f} mol")
    print("-" * 50)


    # --- Step 3: Calculate required moles of n-BuLi ---
    # Equation: Moles n-BuLi = Moles Starting Material * Equivalents
    moles_nBuLi = moles_sm * nBuLi_equivalents
    print("Step 3: Calculate required moles of n-BuLi.")
    print(f"  Equation: Moles_nBuLi = Moles_SM * Equivalents")
    print(f"  Calculation: {moles_sm:.6f} mol * {nBuLi_equivalents:.2f}")
    print(f"  Result: {moles_nBuLi:.6f} mol")
    print("-" * 50)


    # --- Step 4: Calculate the volume of n-BuLi solution to add ---
    # Equation: Volume = Moles / Concentration
    volume_nBuLi_L = moles_nBuLi / nBuLi_concentration_M
    volume_nBuLi_mL = volume_nBuLi_L * 1000
    print("Step 4: Calculate the final volume of n-BuLi solution.")
    print(f"  Equation: Volume = Moles / Concentration")
    print(f"  Calculation: {moles_nBuLi:.6f} mol / {nBuLi_concentration_M:.2f} M")
    print(f"  Result (Liters): {volume_nBuLi_L:.6f} L")
    print(f"  Result (Milliliters): {volume_nBuLi_mL:.3f} mL")
    print("-" * 50)


    print("\n" + "=" * 20 + " FINAL INSTRUCTION " + "=" * 20)
    print(f"For the reaction to be selective, add precisely {volume_nBuLi_mL:.2f} mL of the {nBuLi_concentration_M:.2f} M n-BuLi solution.")
    print("This precision prevents the side reaction that creates the second boron-containing product.")
    print("=" * 61)


if __name__ == '__main__':
    solve_chemistry_problem()