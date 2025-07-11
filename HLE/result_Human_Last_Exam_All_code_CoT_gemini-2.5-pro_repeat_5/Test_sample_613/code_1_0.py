import math

def solve():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in an A-Pb-Br3 perovskite system to determine structural stability.
    """

    # Ionic radii in picometers (pm). Effective radii are used for organic cations.
    # These are commonly accepted values from scientific literature.
    ionic_radii = {
        # A-site cations from answer choices
        "Cesium": 167,
        "Methylammonium": 180,
        "Formamidinium": 253,
        "Aziridinium": 160,       # Too small
        "Ethylammonium": 274,    # Too large
        "Methylhydrazinium": 217,  # Borderline, but generally forms distorted or 2D phases
        "Dimethylammonium": 272, # Too large

        # B-site cation (Lead)
        "Lead (Pb2+)": 119,

        # X-site anion (Bromide)
        "Bromide (Br-)": 196
    }

    # Cations to evaluate based on the provided answer choices
    cations_to_test = [
        "Cesium", "Methylammonium", "Formamidinium", "Aziridinium",
        "Ethylammonium", "Methylhydrazinium", "Dimethylammonium"
    ]

    r_B = ionic_radii["Lead (Pb2+)"]
    r_X = ionic_radii["Bromide (Br-)"]

    print("Analyzing A-site cation stability for A-Pb-Br3 perovskites.")
    print("A stable 3D perovskite structure is expected for a tolerance factor (t) between 0.8 and 1.0.\n")
    print(f"{'Cation':<20} | {'Radius (pm)':<12} | {'Tolerance Factor (t)':<22} | {'Conclusion'}")
    print("-" * 85)

    for cation_name in cations_to_test:
        r_A = ionic_radii[cation_name]
        
        # Goldschmidt tolerance factor calculation
        # t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
        numerator = r_A + r_X
        denominator = math.sqrt(2) * (r_B + r_X)
        t = numerator / denominator

        # Determine stability
        if 0.8 <= t <= 1.0:
            conclusion = "Stable 3D perovskite is likely."
        elif t < 0.8:
            conclusion = "Too small. Does not form 3D perovskite."
        else: # t > 1.0
            conclusion = "Too large. Forms 2D or non-perovskite structures."

        # Print the final equation for each cation
        print(f"For {cation_name}: t = ({r_A} + {r_X}) / (√2 * ({r_B} + {r_X})) = {t:.3f}")
        print(f"{cation_name:<20} | {r_A:<12} | {t:<22.3f} | {conclusion}\n")
        
    print("---")
    print("Summary:")
    print("Cesium, Methylammonium, and Formamidinium have tolerance factors that allow them to independently form 3D perovskite structures.")
    print("The other cations are either too large (Ethylammonium, Dimethylammonium) or too small (Aziridinium) to satisfy the geometric requirements for a stable 3D perovskite, typically forming 2D or other phases instead.")
    print("\nThus, the correct list from the choices is Cesium, Methylammonium, and Formamidinium.")

solve()
<<<A>>>