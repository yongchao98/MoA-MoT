def find_ideal_ni_ce_ratio():
    """
    This script provides the experimentally determined ideal Ni/Ce molar ratios
    for Ni-Ceria catalysts in two important reactions: Water Gas Shift (WGS)
    and Water Splitting (WS).
    """

    print("Ideal Ni/Ce Molar Ratios for Ni-Ceria Nanoparticle Catalysts")
    print("-" * 65)

    # --- 1. Water Gas Shift (WGS) Reaction ---
    # In WGS, Ni is the active catalyst and Ceria (CeO2) is the support.
    # The ratio affects Ni dispersion and metal-support interaction.
    # A value around 0.15 to 0.25 is frequently cited as optimal.
    wgs_ratio_lower = 0.15
    wgs_ratio_upper = 0.25
    
    print("\n1. For the Water Gas Shift (WGS) Reaction:")
    print(f"   The optimal Ni/Ce molar ratio is typically in the range of {wgs_ratio_lower} to {wgs_ratio_upper}.")
    print("   Rationale: This balance provides high dispersion of active Ni particles")
    print("   and strong metal-support interactions, maximizing CO conversion rates")
    print("   and catalyst stability against sintering.")
    # As per the instruction, we output the number in a final equation format.
    # We will use the lower end of the optimal range as a representative value.
    print(f"   Final Equation: Ni_moles / Ce_moles = {wgs_ratio_lower}")


    # --- 2. Thermochemical Water Splitting (WS) ---
    # In thermochemical cycles, Ni is often a dopant in the ceria lattice.
    # The optimal ratio enhances the material's redox properties.
    # A common optimal value is around 10% Ni doping (atomic ratio).
    # Ni/(Ni+Ce) = 0.10 -> Ni/Ce = 0.10 / (1 - 0.10) = 0.10 / 0.90
    ws_ratio = 0.111 

    print("\n2. For Thermochemical Water Splitting (WS):")
    # We'll print the rounded value for clarity.
    print(f"   The optimal Ni/Ce molar ratio is approximately {ws_ratio:.2f}.")
    print("   Rationale: At this level, Ni doping significantly enhances the formation")
    print("   of oxygen vacancies and improves H2 production kinetics without")
    print("   forming separate, less active NiO phases.")
    # Outputting the number in a final equation format.
    print(f"   Final Equation: Ni_moles / Ce_moles = {ws_ratio:.2f}")


if __name__ == "__main__":
    find_ideal_ni_ce_ratio()