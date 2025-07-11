import sys

def calculate_ni_ce_ratio(ni_weight_percent):
    """
    Calculates the Ni/Ce atomic ratio for a Ni/CeO2 catalyst.

    Args:
        ni_weight_percent (float): The weight percent of Nickel in the catalyst.
    """
    # Atomic and Molar Weights (in g/mol)
    ATOMIC_WEIGHT_NI = 58.69
    MOLAR_MASS_CEO2 = 172.12  # (Ce: 140.12 + O2: 2 * 16.00)

    if not 0 < ni_weight_percent < 100:
        print("Error: Nickel weight percent must be between 0 and 100.", file=sys.stderr)
        return

    # Assume a total catalyst mass of 100g for calculation
    mass_ni = float(ni_weight_percent)
    mass_ceo2 = 100.0 - mass_ni

    # Calculate moles of each component
    moles_ni = mass_ni / ATOMIC_WEIGHT_NI
    # Since there is 1 mole of Ce in 1 mole of CeO2, moles of Ce = moles of CeO2
    moles_ce = mass_ceo2 / MOLAR_MASS_CEO2

    # Calculate the final atomic ratio
    atomic_ratio = moles_ni / moles_ce

    print(f"Calculating the Ni/Ce atomic ratio for a catalyst with {ni_weight_percent:.1f} wt% Ni loading.")
    print("-" * 50)
    print(f"Atomic Weight of Ni: {ATOMIC_WEIGHT_NI} g/mol")
    print(f"Molar Mass of CeO2: {MOLAR_MASS_CEO2} g/mol")
    print("\nAssuming 100g of total catalyst:")
    print(f"  Mass of Ni = {mass_ni:.1f} g")
    print(f"  Mass of CeO2 = {mass_ceo2:.1f} g")

    print("\nThe calculation for the atomic ratio is (moles of Ni) / (moles of Ce):")
    # As requested, printing each number in the final equation
    print(f"  Equation: ({mass_ni:.1f} g / {ATOMIC_WEIGHT_NI}) / ({mass_ceo2:.1f} g / {MOLAR_MASS_CEO2})")
    print(f"  Step 1:   ({moles_ni:.4f} moles Ni) / ({moles_ce:.4f} moles Ce)")
    print(f"  Result:   {atomic_ratio:.4f}")
    print("-" * 50)
    print(f"\nThe Ni/Ce atomic ratio is approximately {atomic_ratio:.3f}.")
    print("\nNote: While this is a useful conversion, studies suggest an ideal ratio is often around 0.25 to maximize active sites while preventing nickel particle sintering.")

if __name__ == '__main__':
    # A Ni loading of ~8-10 wt% often corresponds to the cited ideal ratio of ~0.25.
    # We will use 10 wt% as a representative example for this calculation.
    representative_ni_loading = 10.0
    calculate_ni_ce_ratio(representative_ni_loading)